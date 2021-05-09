const express = require('express')
const cors = require('cors')
const app = express()
const port = process.env.PORT || 7777

app.use(cors({ origin: '*' }))

const v1 = require('./v1')
const v2 = require('./v2')
const v3 = require('./v3')

app.get('/', v3)
app.get('/v1', v1)
app.get('/v2', v2)
app.get('/v3', v3)

app.listen(port, () => console.log(`[server] up :${port}`))
