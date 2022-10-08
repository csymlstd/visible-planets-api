
const { Astronomy, Angle, GeographicCoordinates } = require('../astronomy')
const DateTime = require('./utils/datetime')

const NakedEyeObjects = [
    'Sun', 'Moon', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn'
]

const Bodies = Astronomy.Body.filter(body => ['planet', undefined].includes(body.BodyType))

class Sky {

    constructor(opts = {}) {
        const defaults = {
            elevation: 0,
            time: DateTime.create().utc(false).toDate()
        }

        let options = Object.assign({}, defaults, opts)

        this.location = new GeographicCoordinates(options.latitude, options.longitude, options.elevation)
        this.time = Astronomy.DayValue(options.time)
    }

    get(options = {}) {
        let output = []
        Bodies.forEach(body => {
            let item = {
                name: body.Name,
                type: body.BodyType
            }

            const eq = body.EquatorialCoordinates(this.time, this.location)
            const rac = this.getRightAscension(eq.longitude)
            const dec = this.getDeclination(eq.latitude)

            if(item.type) {
                const hc = body.HorizontalCoordinates(this.time, this.location)
                item.aboveHorizon = hc.altitude > 0
            }

            if(options.showCoords) {
                item.rightAscension = rac
                item.declination = dec
            }

            item.nakedEyeObject = NakedEyeObjects.indexOf(body.name) > 0

            delete item.type
            output.push(item)
        })

        if(options.aboveHorizon) {
            output = output.filter(item => item.aboveHorizon === true)
        }

        return output
    }

    getRightAscension(ra) {
        const dms = Angle.DMS(ra)
        if (dms.negative) {
            return null
        }

        const hours   = Number(dms.degrees)
        const minutes = Number(dms.minutes)
        const seconds = Number(dms.seconds)

        return {
            hours, minutes, seconds
        }
    }

    getDeclination(dec) {
        const dms = Angle.DMS(dec)

        const degrees   = Number(dms.degrees)
        const arcminutes = Number(dms.minutes)
        const arcseconds = Number(dms.seconds)

        return {
            degrees: (dms.negative ? degrees * -1 : degrees),
            arcminutes,
            arcseconds
        }
    }

}

module.exports = Sky
