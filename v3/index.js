const Sky = require('./sky');
const DateTime = require('./utils/datetime');
const config = require('../config');

module.exports = (req, res, next) => {
    const latitude = Number(req.query.latitude) || config.defaults.latitude;
    const longitude = Number(req.query.longitude) || config.defaults.longitude;
    const elevation = Number(req.query.elevation) || config.defaults.elevation;
    const time = DateTime.root.utc(req.query.time);
    const aboveHorizon = req.query.aboveHorizon === 'false' ? false : true;
    const showCoords = Boolean(req.query.showCoords);

    const sky = new Sky({
        time: time.toDate(),
        latitude,
        longitude,
        elevation,
    });

    const data = sky.get({ showCoords, aboveHorizon });
    const links = {
        self: 'https://' + req.get('host') + req.originalUrl,
        engine: `https://www.npmjs.com/package/astronomy-engine`,
    };

    return res.json({
        meta: {
            time: time.format(),
            engineVersion: sky.engine_version,
            latitude,
            longitude,
            elevation,
            aboveHorizon
        },
        data,
        links
    });
}
