# Visible Planets API

## Get realtime planetary positions in your sky.
A free JSON API powered by [Don Cross' JS Astronomy Engine](http://cosinekitty.com/astronomy.js).

## Usage

Get a list of planets above the horizon.
```
GET https://visible-planets-api.herokuapp.com?latitude=32&longitude=-98
```

Get a list of planets with their right ascension and declination at the time of the request.
```
GET https://visible-planets-api.herokuapp.com?latitude=32&longitude=-98&showCoords=true
```
