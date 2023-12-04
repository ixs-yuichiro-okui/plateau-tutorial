import { TilesRenderer } from "../src/index.js";
import { GLTFLoader } from "three/examples/jsm/loaders/GLTFLoader.js";
import { DRACOLoader } from "three/examples/jsm/loaders/DRACOLoader.js";

import { CesiumRTCPlugin } from "./CesiumRTCPlugin.js";
import {
	Scene,
	DirectionalLight,
	AmbientLight,
	WebGLRenderer,
	PerspectiveCamera,
	Box3,
	OrthographicCamera,
	sRGBEncoding,
	Group,
	ShaderMaterial,
	MeshStandardMaterial,
	PCFSoftShadowMap,
	Sphere,
	GridHelper,
	AxesHelper,
	MeshLambertMaterial,
	CircleGeometry,
	MeshBasicMaterial,
	Mesh,
	Vector2,
	Vector3,
	Matrix3,
	Matrix4,
	Quaternion,
	Raycaster,
	Plane,
	BoxGeometry,
	PlaneGeometry,
	CameraHelper,
	Box3Helper,
	ArrowHelper,
	EdgesGeometry,
	LineSegments,
	WireframeGeometry,
	BufferGeometry,
	BufferAttribute,
	Vector4
} from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { VertexNormalsHelper } from "three/examples/jsm/helpers/VertexNormalsHelper.js";
import { GeodeticTransformer } from "./geodetic-transform.js";

import { GUI } from "three/examples/jsm/libs/lil-gui.module.min.js";
import Stats from "three/examples/jsm/libs/stats.module.js";
import EcefTransform from "./EcefTransform.js";
import { Ecef2Wgs84, euclideanTransform3D } from "./EcefTransform.js";

let camera, controls, scene, renderer, tiles, orthoCamera;
let prefetchCamera;
let offsetParent, box, matrix, sphereCentroid, dirLight, statsContainer;
let raycaster;
let groundParent;
let stats;
let interactiveModel = [],
	buildingPickList = [];
let masterTileset;
let testHelperPanel;
let focusAxes;

const DEFAULT = 0;
const GRADIENT = 1;
const TOPOGRAPHIC_LINES = 2;
const LIGHTING = 3;
const params = {
	"material": DEFAULT,
	"orthographic": false,
	"picking_building": false,
	"rebuild": initTiles
};

const gradientShader = {
	vertexShader: /* glsl */ `
		varying vec3 wPosition;
		void main() {

			#include <begin_vertex>
			#include <project_vertex>
			wPosition = ( modelMatrix * vec4( transformed, 1.0 ) ).xyz;

		}
	`,
	fragmentShader: /* glsl */ `
		varying vec3 wPosition;
		void main() {

			float minVal = - 30.0;
			float maxVal = 30.0;

			float val = ( wPosition.y - minVal ) / ( maxVal - minVal );

			vec4 color1 = vec4( 0.149, 0.196, 0.219, 1.0 ) * 0.5;
			vec4 color2 = vec4( 1.0 );

			gl_FragColor = mix( color1, color2, val );

		}
	`
};

const topoShader = {
	extensions: {
		derivatives: true
	},
	vertexShader: /* glsl */ `
		varying vec3 wPosition;
		varying vec3 vViewPosition;
		void main() {

			#include <begin_vertex>
			#include <project_vertex>
			wPosition = ( modelMatrix * vec4( transformed, 1.0 ) ).xyz;
			vViewPosition = - mvPosition.xyz;

		}
	`,
	fragmentShader: /* glsl */ `
		varying vec3 wPosition;
		varying vec3 vViewPosition;
		void main() {

			// lighting
			vec3 fdx = vec3( dFdx( wPosition.x ), dFdx( wPosition.y ), dFdx( wPosition.z ) );
			vec3 fdy = vec3( dFdy( wPosition.x ), dFdy( wPosition.y ), dFdy( wPosition.z ) );
			vec3 worldNormal = normalize( cross( fdx, fdy ) );

			float lighting =
				0.4 +
				clamp( dot( worldNormal, vec3( 1.0, 1.0, 1.0 ) ), 0.0, 1.0 ) * 0.5 +
				clamp( dot( worldNormal, vec3( - 1.0, 1.0, - 1.0 ) ), 0.0, 1.0 ) * 0.3;

			// thickness scale
			float upwardness = dot( worldNormal, vec3( 0.0, 1.0, 0.0 ) );
			float yInv = clamp( 1.0 - abs( upwardness ), 0.0, 1.0 );
			float thicknessScale = pow( yInv, 0.4 );
			thicknessScale *= 0.25 + 0.5 * ( vViewPosition.z + 1.0 ) / 2.0;

			// thickness
			float thickness = 0.01 * thicknessScale;
			float thickness2 = thickness / 2.0;
			float m = mod( wPosition.y, 3.0 );

			// soften edge
			float center = thickness2;
			float dist = clamp( abs( m - thickness2 ) / thickness2, 0.0, 1.0 );

			vec4 topoColor = vec4( 0.149, 0.196, 0.219, 1.0 ) * 0.5;
			gl_FragColor = mix( topoColor * lighting, vec4( lighting ), dist );

		}
	`
};

init();
animate();

window.calculateOffset = () => {
	const list = [];
	scene.traverse((c) => {
		if (c.isMesh) {
			list.push(c.getWorldPosition(new Vector3()));
		}
	});

	return list;
};

function updateMaterial(scene) {
	const materialIndex = parseFloat(params.material);
	scene.traverse((c) => {
		if (c.isMesh) {
			c.material.dispose();

			switch (materialIndex) {
				case DEFAULT:
					c.material = c.originalMaterial;
					c.material.side = 2;
					c.receiveShadow = false;
					c.castShadow = false;
					break;
				case GRADIENT:
					c.material = new ShaderMaterial(gradientShader);
					c.material.side = 2;
					c.receiveShadow = false;
					c.castShadow = false;
					break;
				case TOPOGRAPHIC_LINES:
					c.material = new ShaderMaterial(topoShader);
					c.material.side = 2;
					c.material.flatShading = true;
					c.receiveShadow = false;
					c.castShadow = false;
					break;
				case LIGHTING:
					c.material = new MeshStandardMaterial();
					c.material.side = 2;
					c.receiveShadow = true;
					c.castShadow = true;
			}
		}
	});
}

function convertWgs84ToLatLng(wgs84) {
	const [x, y] = wgs84;
	const lng = x * (180 / Math.PI);
	const lat = (Math.atan(Math.exp(y)) * 360) / Math.PI - 90;
	return [lat, lng];
}

function createGeometryByBatchId(geometry, batchID) {
	const rawBatchID = geometry.getAttribute("_batchid").array;

	const newIndices = [];
	for (let i = 0; i < rawBatchID.length; i++) {
		if (rawBatchID[i] === batchID) {
			newIndices.push(i);
		}
	}

	const rawPosition = geometry.getAttribute("position").array;
	const newPosition = [];
	for (let i = 0; i < newIndices.length; i++) {
		const ptr = newIndices[i] * 3;

		newPosition.push(rawPosition[ptr], rawPosition[ptr + 1], rawPosition[ptr + 2]);
		newIndices[i] = i;
	}

	const newGeometry = new BufferGeometry();
	newGeometry.setAttribute("position", new BufferAttribute(new Float32Array(newPosition), 3));
	newGeometry.setIndex(new BufferAttribute(new Uint16Array(newIndices), 1));
	newGeometry.computeBoundingBox();
	newGeometry.computeBoundingSphere();

	const newMaterial = new MeshBasicMaterial({ color: 0x00ff00 });
	const newMesh = new Mesh(newGeometry, newMaterial);

	return newMesh;
}

function onLoadModel(model) {
	const mesh = model.children.filter((c) => c.isMesh)[0];

	if (mesh) {
		const newMesh = createGeometryByBatchId(mesh.geometry, 56);
		if (newMesh) {
			newMesh.roi = true;
			model.add(newMesh);

			const boundingBox = newMesh.geometry.boundingBox;
			const boundingBoxHelper = new Box3Helper(boundingBox, 0xffff00);
			model.add(boundingBoxHelper);

			model.matrixWorldNeedsUpdate = true;
		}
	}

	model.traverse((c) => {
		if (c.isMesh) {
			c.originalMaterial = c.material;
		}
	});

	updateMaterial(model);
}

function onDisposeModel(scene) {
	scene.traverse((c) => {
		if (c.isMesh) {
			c.material.dispose();
		}
	});

	interactiveModel.length = 0;
}

function prefetchTiles() {}

function initTiles() {
	if (tiles) {
		tiles.group.parent.remove(tiles.group);
		tiles.dispose();
	}

	prefetchCamera = new OrthographicCamera();

	//const url = "https://plateau.geospatial.jp/main/data/3d-tiles/bldg/14130_kawasaki/notexture/tileset.json";
	const url = "	https://plateau.geospatial.jp/main/data/3d-tiles/bldg/13100_tokyo/13101_chiyoda-ku/texture/tileset.json";
	//const url = "https://plateau.geospatial.jp/main/data/3d-tiles/bldg/14100_yokohama/notexture/tileset.json";
	//const url = "	https://plateau.geospatial.jp/main/data/3d-tiles/bldg/14150_sagamihara/notexture/tileset.json";
	//const url = "https://plateau.geospatial.jp/main/data/3d-tiles/bldg/14201_yokosuka/notexture/tileset.json";

	tiles = new TilesRenderer(url);
	tiles.optimizeRaycast = false;
	tiles.setCamera(prefetchCamera);
	tiles.setResolutionFromRenderer(prefetchCamera, renderer);

	tiles.onLoadTileSet = (tileset) => {
		tiles.maxDepth = 2;

		const sphere = tileset.root.cached.sphere;
		const dir = sphere.center.clone().normalize().multiplyScalar(-1);
		const quaternion = new Quaternion();
		quaternion.setFromUnitVectors(new Vector3(0, 0, -1), dir);

		prefetchCamera.left = -sphere.radius;
		prefetchCamera.right = sphere.radius;
		prefetchCamera.bottom = -sphere.radius;
		prefetchCamera.top = sphere.radius;
		prefetchCamera.near = 0.1;
		prefetchCamera.far = sphere.center.distanceTo(new Vector3()) + sphere.radius;
		prefetchCamera.position.set(sphere.center.x, sphere.center.y, sphere.center.z).add(sphere.radius, sphere.radius, sphere.radius);
		prefetchCamera.lookAt(new Vector3(0, 0, 0));

		prefetchCamera.updateMatrixWorld(true);
		prefetchCamera.updateProjectionMatrix();

		masterTileset = tileset;
		//tiles.group.rotateZ(Math.PI);
	};

	const dracoLoader = new DRACOLoader();
	dracoLoader.setDecoderConfig({ type: "js" });
	dracoLoader.setDecoderPath("https://www.gstatic.com/draco/v1/decoders/");

	const gltfLoader = new GLTFLoader(tiles.manager);
	tiles.manager.addHandler(/\.gltf$/, gltfLoader);
	gltfLoader.setDRACOLoader(dracoLoader);
	gltfLoader.register((parser) => new CesiumRTCPlugin(parser));

	groundParent = new Group();

	///////////////////tiles.group to JGD///////////////////////////////
	var cnt = 0;
	function tiles_group_to_JGD(model) {
		const children = masterTileset.root.children;
		if (!children.every((element) => element.cached.scene)) {
			return;
		}

		const filterChildren = children.reduce((acc, element) => {
			acc.push(element.cached.scene);
			return acc;
		}, []);

		if (tiles.transform) {
			return;
		}

		// change camera
		tiles.deleteCamera(prefetchCamera);
		tiles.setCamera(camera);
		tiles.setResolutionFromRenderer(camera, renderer);
		tiles.update();

		const createOriginParam = (mesh) => {
			const subMesh = createGeometryByBatchId(mesh.geometry, 0);
			mesh.add(subMesh);

			const subMeshCenter = subMesh.geometry.boundingSphere.center.clone().applyMatrix4(mesh.parent.matrix.clone());

			return {
				position0: mesh.position,
				position: subMeshCenter,
				subMeshCenter: subMeshCenter
			};
		};

		//const c1 = tiles.group.children[0];

		const c1 = filterChildren[0];
		const test_c1 = filterChildren[0];
		const c1Param = c1; //createOriginParam(c1.children[0]);
		//const c1Geometry = createGeometryByBatchId(c1.children[0].geometry, 0);

		//const c2 = tiles.group.children[1];
		const c2 = filterChildren[1];
		const c2Param = c2; //createOriginParam(c2.children[0]);
		//const c2Geometry = createGeometryByBatchId(c2.children[0].geometry, 0);

		//const c3 = tiles.group.children[2];
		const c3 = filterChildren[2];
		const c3Param = c3; //createOriginParam(c3.children[0]);
		//const c3Geometry = createGeometryByBatchId(c3.children[0].geometry, 0);

		/*
		const c1 = filterChildren[0];
		const test_c1 = filterChildren[0];
		const c1Param = createOriginParam(c1.children[0]);

		const c2 = filterChildren[1];
		const c2Param = createOriginParam(c2.children[0]);

		const c3 = filterChildren[2];
		const c3Param = createOriginParam(c3.children[0]);
		*/
		console.log(c1, c1Param);
		console.log(c2, c2Param);
		console.log(c3, c3Param);

		const CRTC2JGD = (position) => {
			const obj = Ecef2Wgs84(position);
			const input = {};
			input.longitudeLatitude = {};
			input.longitudeLatitude.x = obj.longitude;
			input.longitudeLatitude.y = obj.latitude;

			const conv = GeodeticTransformer.WorldGeodetic2JapaneseDatum(input, 9);
			conv.z = obj.h;
			return new Vector3(conv.xy.x, conv.xy.y, obj.h);
		};

		const veclistToMatrix = (vecList) => {
			const result = [];
			for (let vec of vecList) {
				result.push([vec.x, vec.y, vec.z]);
			}
			return result;
		};

		console.log(c1Param.position, CRTC2JGD(c1Param.position));
		console.log(c2Param.position, CRTC2JGD(c2Param.position));
		console.log(c3Param.position, CRTC2JGD(c3Param.position));

		const A = veclistToMatrix([c1Param.position, c2Param.position, c3Param.position]);
		const B = veclistToMatrix([CRTC2JGD(c1Param.position), CRTC2JGD(c2Param.position), CRTC2JGD(c3Param.position)]);

		const Rt = euclideanTransform3D(A, B);
		const T = new Matrix4();
		T.set(Rt.R[0][0], Rt.R[0][1], Rt.R[0][2], Rt.t[0], Rt.R[1][0], Rt.R[1][1], Rt.R[1][2], Rt.t[1], Rt.R[2][0], Rt.R[2][1], Rt.R[2][2], Rt.t[2], 0, 0, 0, 1);

		tiles.group.applyMatrix4(T);

		const R = new Matrix3();
		R.set(Rt.R[0][0], Rt.R[0][1], Rt.R[0][2], Rt.R[1][0], Rt.R[1][1], Rt.R[1][2], Rt.R[2][0], Rt.R[2][1], Rt.R[2][2]);
		const t = new Vector3(Rt.t[0], Rt.t[1], Rt.t[2]);

		console.log(test_c1.position);
		console.log(R);
		console.log(test_c1.position.clone().applyMatrix3(R.transpose()));
		console.log(t);
		console.log(test_c1.position.clone().applyMatrix3(R.transpose()).add(t));

		tiles.transform = T;

		if (tiles.getOrientedBounds(box, matrix)) {
			box.min.z = box.max.z = Math.min(box.min.z, box.max.z);

			const center = new Vector3();
			box.applyMatrix4(matrix);
			box.getCenter(center);

			focusTiles(center, matrix);
		}
		tiles.maxDepth = Infinity;

		if (tiles.getBoundingSphere(sphereCentroid)) {
			//focusTiles(sphereCentroid.center);
		}
	}

	////////////////////////////////////////////////////////////////////

	tiles.onLoadModel = (model) => {
		onLoadModel(model);
		tiles_group_to_JGD(model);
	};
	tiles.onDisposeModel = onDisposeModel;

	scene.add(tiles.group);

	window.tiles = tiles;
	window.tilesObject = tiles.group;
}

function focusTiles(center) {
	const centralPosition = center.clone().applyMatrix4(tiles.group.matrix);
	const centralQuaternion = new Quaternion();
	centralQuaternion.setFromRotationMatrix(tiles.group.matrix);

	const gridHelper = new GridHelper(100000, 1000);
	tiles.group.add(gridHelper);

	controls.target.copy(centralPosition);
	camera.quaternion.copy(centralQuaternion);

	if (!focusAxes) {
		focusAxes = new AxesHelper(2000);
		focusAxes.position.copy(centralPosition);
		scene.add(focusAxes);
	}
}

function init() {
	scene = new Scene();

	// primary camera view
	renderer = new WebGLRenderer({ antialias: true });
	renderer.setPixelRatio(window.devicePixelRatio);
	renderer.setSize(window.innerWidth, window.innerHeight);
	renderer.setClearColor(0x151c1f);
	renderer.shadowMap.enabled = true;
	renderer.shadowMap.type = PCFSoftShadowMap;
	renderer.outputEncoding = sRGBEncoding;

	document.body.appendChild(renderer.domElement);

	camera = new PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.1, 10000000);

	orthoCamera = new OrthographicCamera();

	// controls
	controls = new OrbitControls(camera, renderer.domElement);
	controls.screenSpacePanning = false;
	controls.minDistance = 1;
	controls.maxDistance = 20000;

	// lights
	dirLight = new DirectionalLight(0xffffff, 1.25);
	dirLight.position.set(1, 2, 3).multiplyScalar(40);
	dirLight.castShadow = true;
	dirLight.shadow.bias = -0.01;
	dirLight.shadow.mapSize.setScalar(2048);

	// raycaster
	raycaster = new Raycaster();
	raycaster.firstHitOnly = true;
	let selectModel = { id: -1 };
	const selectMaterial = new MeshLambertMaterial({
		transparent: true,
		opacity: 0.6,
		color: 0xff00ff,
		depthTest: false
	});

	renderer.domElement.ondblclick = (event) => translate(event, selectMaterial, selectModel, camera);

	const shadowCam = dirLight.shadow.camera;
	shadowCam.left = -200;
	shadowCam.bottom = -200;
	shadowCam.right = 200;
	shadowCam.top = 200;
	shadowCam.updateProjectionMatrix();

	const axesHelper = new AxesHelper(100000);
	scene.add(axesHelper);

	scene.add(dirLight);

	const ambLight = new AmbientLight(0xffffff, 0.05);
	scene.add(ambLight);

	box = new Box3();
	matrix = new Matrix4();
	sphereCentroid = new Sphere();

	// offsetParent = new Group();
	// scene.add(offsetParent);

	window.scene = scene;
	window.camera = camera;
	window.controls = controls;
	window.offsetParent = offsetParent;

	initTiles();

	onWindowResize();
	window.addEventListener("resize", onWindowResize, false);

	// GUI
	const gui = new GUI();
	gui.width = 300;
	gui.add(params, "orthographic");
	gui.add(params, "picking_building");
	gui.add(params, "material", { DEFAULT, GRADIENT, TOPOGRAPHIC_LINES, LIGHTING }).onChange(() => {
		tiles.forEachLoadedModel(updateMaterial);
	});
	gui.add(params, "rebuild");
	gui.open();

	// Stats
	stats = new Stats();
	stats.showPanel(0);
	document.body.appendChild(stats.dom);

	statsContainer = document.createElement("div");
	statsContainer.style.position = "absolute";
	statsContainer.style.top = 0;
	statsContainer.style.left = 0;
	statsContainer.style.color = "white";
	statsContainer.style.width = "100%";
	statsContainer.style.textAlign = "center";
	statsContainer.style.padding = "5px";
	statsContainer.style.pointerEvents = "none";
	statsContainer.style.lineHeight = "1.5em";
	document.body.appendChild(statsContainer);
}

function calculateMousePosition(event) {
	const mouse = new Vector2();

	// スクリーン上のマウスの位置を計算する
	const bounds = renderer.domElement.getBoundingClientRect();

	const x1 = event.clientX - bounds.left;
	const x2 = bounds.right - bounds.left;
	mouse.x = (x1 / x2) * 2 - 1;

	const y1 = event.clientY - bounds.top;
	const y2 = bounds.bottom - bounds.top;
	mouse.y = -(y1 / y2) * 2 + 1;

	return mouse;
}

function cast(position, camera) {
	raycaster.setFromCamera(position, camera);

	return raycaster.intersectObject(tiles.group, true);
}

function translate(event, material, model, camera = camera, is_debug = false) {
	const mouse = calculateMousePosition(event);
	const found = cast(mouse, camera);

	if (found && found.length > 0) {
		const roiList = found.filter((c) => c.object.roi);

		if (roiList.length == 0) {
			return;
		}

		const roiObject = roiList[0].object;

		console.log(roiObject);

		const meshMatrix = groundParent.matrix.clone().multiply(tiles.group.matrix).multiply(roiObject.parent.matrix);
		const meshRotation = new Quaternion();
		meshRotation.setFromRotationMatrix(meshMatrix);

		const origin = new Vector3().applyMatrix4(meshMatrix).multiplyScalar(-1);
		offsetParent.position.copy(origin);

		const geometryPosition = new Vector3();
		geometryPosition.copy(roiObject.geometry.boundingSphere.center.clone());
		geometryPosition.applyQuaternion(meshRotation);

		const diffDir = geometryPosition.clone().normalize();
		const arrowHelper = new ArrowHelper(diffDir, origin.multiplyScalar(-1), geometryPosition.distanceTo(new Vector3()));

		offsetParent.add(arrowHelper);
		offsetParent.position.add(geometryPosition.multiplyScalar(-1));
	}
}

function onWindowResize() {
	camera.aspect = window.innerWidth / window.innerHeight;
	renderer.setPixelRatio(window.devicePixelRatio);
	renderer.setSize(window.innerWidth, window.innerHeight);
	camera.updateProjectionMatrix();

	updateOrthoCamera();
}

function updateOrthoCamera() {
	orthoCamera.position.copy(camera.position);
	orthoCamera.rotation.copy(camera.rotation);

	const scale = camera.position.distanceTo(controls.target) / 2.0;
	const aspect = window.innerWidth / window.innerHeight;
	orthoCamera.left = -aspect * scale;
	orthoCamera.right = aspect * scale;
	orthoCamera.bottom = -scale;
	orthoCamera.top = scale;
	orthoCamera.near = camera.near;
	orthoCamera.far = camera.far;
	orthoCamera.updateProjectionMatrix();
}

let animateFlag = false;

function animate() {
	requestAnimationFrame(animate);

	/*
	if (params.orthographic) {
		tiles.deleteCamera(camera);
		tiles.setCamera(orthoCamera);
		tiles.setResolutionFromRenderer(orthoCamera, renderer);
	} else {
		tiles.deleteCamera(orthoCamera);
		tiles.setCamera(camera);
		tiles.setResolutionFromRenderer(camera, renderer);
	}
	*/

	/*
	groundParent.rotation.set(0, 0, 0);
	if (params.up === "-Z") {
		groundParent.rotation.x = Math.PI / 2;
	}
	groundParent.updateMatrixWorld(true);
	*/

	const correctInclination = (center) => {
		if (!animateFlag) {
			camera.position.copy(center);

			//controls.target.copy(center);
			//controls.position0.copy(center);

			// tiles.group.position.copy(center);
			// tiles.group.position.multiplyScalar(-1);
			// tiles.group.updateMatrixWorld();

			animateFlag = true;
		}

		// const obj = Ecef2Wgs84(center);
		// const input = {};
		// input.longitudeLatitude = {};
		// input.longitudeLatitude.x = obj.longitude;
		// input.longitudeLatitude.y = obj.latitude;
		//console.log(GeodeticTransformer.WorldGeodetic2JapaneseDatum(input, 9));

		/*
		const direction = center.clone().normalize();
		const up = new Vector3(0, 1, 0);
		groundParent.quaternion.setFromUnitVectors(direction, up);
		groundParent.updateMatrixWorld();
			*/
	};

	if (tiles.getOrientedBounds(box, matrix)) {
		box.min.z = box.max.z = Math.min(box.min.z, box.max.z);

		const center = new Vector3();
		box.applyMatrix4(matrix);
		box.getCenter(center);

		correctInclination(center);
	}

	if (tiles.getBoundingSphere(sphereCentroid)) {
		//correctInclination(sphereCentroid.center);
		//console.log(Ecef2Wgs84(sphereCentroid.center));
	}

	// update tiles
	window.tiles = tiles;
	camera.updateMatrixWorld();
	orthoCamera.updateMatrixWorld();
	//prefetchCamera.updateMatrixWorld();
	tiles.update();

	render();
	stats.update();
}

function render() {
	updateOrthoCamera();

	statsContainer.innerText = `Geometries: ${renderer.info.memory.geometries} ` + `Textures: ${renderer.info.memory.textures} ` + `Programs: ${renderer.info.programs.length} `;

	//renderer.render(scene, prefetchCamera);
	renderer.render(scene, params.orthographic ? orthoCamera : camera);
}
