const Sky = require('./sky')
const DateTime = require('./utils/datetime')

module.exports = (req, res, next) => {
    const latitude = Number(req.query.latitude) || 32.78
    const longitude = Number(req.query.longitude) || -96.84
    const elevation = Number(req.query.elevation) || 0
    const time = DateTime.root.utc(req.query.time)
    const aboveHorizon = req.query.aboveHorizon === 'false' ? false : true
    const showCoords = Boolean(req.query.showCoords)

    const sky = new Sky({
        time: time.toDate(),
        latitude,
        longitude,
        elevation,
    })

    const data = sky.get({ showCoords, aboveHorizon })
    const links = {
        self: req.protocol + '://' + req.get('host') + req.originalUrl,
        engine: `https://www.npmjs.com/package/astronomy-engine`,
    }

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
    })
}
