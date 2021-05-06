const Sky = require('./sky')

module.exports = (req, res, next) => {
    const latitude = Number(req.query.latitude) || 32.78
    const longitude = Number(req.query.longitude) || -96.84
    const elevation = Number(req.query.elevation) || null

    const sky = new Sky({
        latitude, longitude, elevation
    })

    return res.json(sky.get({ showCoords: req.query.showCoords }))
}