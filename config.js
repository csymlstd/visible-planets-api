const fs = require('fs');
const path = require('path');

const getAstronomyEngineVersion = () => {
    try {
        const enginePath = path.resolve('./node_modules/astronomy-engine/package.json');
        const info = JSON.parse(fs.readFileSync(enginePath, {encoding: 'utf8'}));
        return info.version;
    } catch(err) {
        return 'unknown';
    }
}

module.exports = {
    engineVersion: getAstronomyEngineVersion(),
    defaults: {
        latitude: 28.627222,
        longitude: -80.620833,
        elevation: 0,
    }
}
