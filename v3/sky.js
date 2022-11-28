const { AstroTime, Body, Angle, Observer, Equator, Horizon, Illumination, MoonPhase, Constellation, Elongation } = require('astronomy-engine');
const DateTime = require('./utils/datetime');
const config = require('../config');

const NakedEyeObjects = [
    Body.Sun, Body.Moon, Body.Mercury, Body.Venus, Body.Earth, Body.Mars, Body.Jupiter, Body.Saturn,
];

const IncludedBodies = [
    Body.Sun, Body.Moon, Body.Mercury, Body.Venus, Body.Mars,
    Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto
];

class Sky {

    constructor(opts = {}) {
        const defaults = {
            elevation: 0,
            time: DateTime.create().utc(false).toDate()
        }

        let options = Object.assign({}, defaults, opts)
        this.engine_version = config.engineVersion;
        this.location = new Observer(options.latitude, options.longitude, options.elevation);
        this.time = new AstroTime(options.time);
    }

    getDMS(x) {
        var a = {};

        a.negative = (x < 0);
        if (a.negative) {
            x = -x;
        }

        a.degrees = Math.floor(x);
        x = 60.0 * (x - a.degrees);
        a.arcminutes = Math.floor(x);
        x = 60.0 * (x - a.arcminutes);
        a.arcseconds = Math.round(10.0 * x) / 10.0;   // Round to the nearest tenth of an arcsecond.

        return a;
    }

    getHMS(x) {
        var a = {};

        a.negative = (x < 0);
        if (a.negative) {
            x = -x;
        }

        a.hours = Math.floor(x);
        x = 60.0 * (x - a.hours);
        a.minutes = Math.floor(x);
        x = 60.0 * (x - a.minutes);
        a.seconds = Math.round(10.0 * x) / 10.0;   // Round to the nearest tenth of an arcsecond.

        return a;
    }

    get(options = {}) {
        let output = []
        Object.keys(Body).filter(key => IncludedBodies.includes(key)).forEach(body => {
            let item = {
                name: body,
            }

            const eq = new Equator(Body[body], this.time, this.location, true, false)
            const hc = new Horizon(this.time, this.location, eq.ra, eq.dec, 'normal')
            const constellation = new Constellation(eq.ra, eq.dec)

            item.constellation = constellation.name
            item.rightAscension = this.getHMS(hc.ra)
            item.declination = this.getDMS(hc.dec)
            item.rightAscension.raw = hc.ra
            item.declination.raw = hc.dec
            item.altitude = hc.altitude
            item.azimuth = hc.azimuth
            item.aboveHorizon = hc.altitude > 0

            if(body == Body.Moon) {
                item.phase = (MoonPhase(this.time) / 360) * 100
            }

            const illumination = new Illumination(Body[body], this.time)
            item.magnitude = illumination.mag

            item.nakedEyeObject = NakedEyeObjects.indexOf(body) > -1

            if((item.aboveHorizon && options.aboveHorizon) || !options.aboveHorizon) {
                output.push(item)
            }
        })

        return output
    }

}

module.exports = Sky
