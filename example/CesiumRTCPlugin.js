export class CesiumRTCPlugin {
	constructor(parser) {
		this.parser = parser;
	}

	afterRoot(result) {
		if (this.parser.json.extensions?.CESIUM_RTC?.center != null) {
			const center = this.parser.json.extensions.CESIUM_RTC.center;
			result.scene.position.set(...center);
		}
	}
}
