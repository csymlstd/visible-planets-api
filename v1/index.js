const Sky = require('./sky');
const config = require('../config');

module.exports = (req, res, next) => {
    const latitude = Number(req.query.latitude) || config.defaults.latitude;
    const longitude = Number(req.query.longitude) || config.defaults.longitude;
    const elevation = Number(req.query.elevation) || config.defaults.elevation;

    const sky = new Sky({
        latitude, longitude, elevation
    });

    return res.json(sky.get({ showCoords: req.query.showCoords }));
}
