# Visible Planets API

## Get realtime planetary positions in your sky.
A free JSON API powered by [Don Cross' JS Astronomy Engine](http://cosinekitty.com/astronomy.js).

## Usage

Get a list of planets (and Moon) above the horizon.
```
GET https://visible-planets-api.herokuapp.com/v2?latitude=32&longitude=-98
```

Get a list of planets (and Moon) with their declination and right ascension coordinates.
```
GET https://visible-planets-api.herokuapp.com/v2?latitude=32&longitude=-98&showCoords=true
```

## Query Parameters

| Param | Default Value | Description |
| ----- | ------------- | ----------- |
| latitude | 32.78 | Latitude of observer |
| longitude | -96.84 | Longitude of observer |
| elevation | 0 | Elevation of observer in meters above sea level |
| time | null | Time of observation in [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601) format, defaults to time of request |
| showCoords | false | Display declination and right ascension of each body, expects true or false

## Changelog

### v2 - 2021-05-04
- Changed declination response properties from hours to degrees, minutes to arcminutes, seconds to arcseconds,
- Changed declination and right ascension response values from strings to numbers
- Added query parameter to set time of observation using ISO 8601 format
- Fixed typo in `rightAscension` response parameter
- Response now follows JSON:API spec and includes parameters used to generate sky