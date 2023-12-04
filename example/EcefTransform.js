const math = require("mathjs");
import { SVD } from "svd-js";

const ellipse_a = 6378137.0;
const flattening = 1 / 298.257223563;
const ellipse_b = ellipse_a * (1 - flattening);
const ecef_e2 = 2 * flattening - flattening ** 2;
const ecef_ed2 = (ecef_e2 * ellipse_a ** 2) / ellipse_b ** 2;
const coefRad2Deg = 180 / Math.PI;

export function Ecef2Wgs84(ecef) {
	const p = Math.sqrt(ecef.x ** 2 + ecef.y ** 2);
	const theta = Math.atan2(ecef.z * ellipse_a, p * ellipse_b);

	const latitude = Math.atan2(ecef.z + ecef_ed2 * ellipse_b * Math.sin(theta) ** 3, p - ecef_e2 * ellipse_a * Math.cos(theta) ** 3);
	const longitude = Math.atan2(ecef.y, ecef.x);

	const N = (latitude) => {
		return ellipse_a / Math.sqrt(1 - ecef_e2 * Math.sin(latitude) ** 2);
	};

	const h = p / Math.cos(latitude) - N(latitude);

	return {
		latitude: latitude * coefRad2Deg,
		longitude: longitude * coefRad2Deg,
		h: h
	};
}

export function euclideanTransform3D(A, B, args = {}) {
	const default_config = {
		eps: args.eps ?? 1e-10,
		logging: args.logging ?? false
	};

	// A, B are arrays of arrays, e.g. [[x1, y1, z1], [x2, y2, z2], ...]
	if (A.length !== B.length) {
		throw new Error("A and B must have the same length.");
	}

	const N = A.length;

	const centroidA = math.mean(A, 0);

	const centroidB = math.mean(B, 0);

	const AA = A.map((point) => math.subtract(point, centroidA));

	const BB = B.map((point) => math.subtract(point, centroidB));

	const H = math.multiply(math.transpose(AA), BB);

	if (default_config.logging) console.log(H);
	const svdResult = SVD(H);

	const S = svdResult.q;
	const checkDim = S.filter((value) => Math.abs(value) < default_config.eps);
	if (checkDim.length >= 2) {
		if (default_config.logging) console.warn("Eigenvalue with SVD was insufficient.");
	}

	const U = svdResult.u;
	const Vt = math.transpose(svdResult.v);

	let R = math.multiply(math.transpose(Vt), math.transpose(U));

	const testMat = [
		[0.4131232159761867, -0.3499306784349872, 0.8407602088060079],
		[-0.6203769102569903, -0.7840097134441982, -0.02147692820531234],
		[0.6666796064388673, -0.5127156029566249, -0.5409815272662658]
	];
	console.log(testMat);
	const dettest = math.det(testMat);
	console.log(testMat);

	if (math.det(R) < 0) {
		if (default_config.logging) console.log("sign");

		Vt[2] = math.multiply(Vt[2], -1);

		R = math.multiply(math.transpose(Vt), math.transpose(U));

		if (default_config.logging) console.log("new R", R);
	}

	const t = math.subtract(centroidB, math.multiply(R, centroidA));

	return { R, t };
}
