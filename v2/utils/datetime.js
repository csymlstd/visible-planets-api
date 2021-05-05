const dayjs = require('dayjs')
const utc = require('dayjs/plugin/utc')
const timezone = require('dayjs/plugin/timezone')
const customParseFormat = require('dayjs/plugin/customParseFormat')
const advancedFormat = require('dayjs/plugin/advancedFormat')

dayjs.extend(utc)
dayjs.extend(timezone)
dayjs.extend(customParseFormat)
dayjs.extend(advancedFormat)

module.exports = class DateTime {

    static create(date, formatStr) {
        return dayjs(date, formatStr)
    }

    static get root() {
        return dayjs
    }

}