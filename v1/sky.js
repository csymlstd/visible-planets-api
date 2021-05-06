
const { Astronomy, Angle, GeographicCoordinates } = require('../astronomy')

const NakedEyeObjects = [
    'Sun', 'Moon', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn'
]

const Bodies = Astronomy.Body.filter(body => ['planet', undefined].includes(body.BodyType))

class Sky {
    
    constructor(opts = {}) {
        const defaults = {
            elevation: null,
            time: new Date()
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
                item.rightAscenson = rac
                item.declination = dec
            }
            
            item.nakedEyeObject = NakedEyeObjects.indexOf(body.name) > 0
            
            delete item.type
            output.push(item)
        })

        output = output.filter(item => item.aboveHorizon === true)

        return output
    }

    getRightAscension(ra) {
        const dms = Angle.DMS(ra)
        if (dms.negative) {
            return null
        }
        
        const hours   = dms.degrees.toFixed()
        const minutes = dms.minutes.toFixed()
        const seconds = dms.seconds.toFixed(1)

        return {
            hours, minutes, seconds
        }
    }

    getDeclination(dec) {
        const dms = Angle.DMS(dec)

        const hours   = dms.degrees.toFixed(1)
        const minutes = dms.minutes.toFixed(1)
        const seconds = dms.seconds.toFixed(1)

        return {
            hours: (dms.negative ? hours * -1 : hours).toString(),
            minutes,
            seconds
        }
    }

}

module.exports = Sky