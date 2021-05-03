/*
    astronomy.js  -  by Don Cross  -  http://cosinekitty.com
    
    General-purpose astronomy calculation routines.
    Designed for use by html, svg, and cscript/wscript.
*/

//----------------------------------------------------------------------------------------------
// Numeric constants...

var METERS_PER_ASTRONOMICAL_UNIT = 1.4959787e+11;
var METERS_PER_EARTH_EQUATORIAL_RADIUS = 6378140.0;
var EARTH_RADII_PER_ASTRONOMICAL_UNIT = METERS_PER_ASTRONOMICAL_UNIT / METERS_PER_EARTH_EQUATORIAL_RADIUS;      // 23454.78

//----------------------------------------------------------------------------------------------
// Custom trigonometry and math functions.
// We do most of our math in angular degrees and sideral hours, instead of radians.

function AngleClass()
{
    this.DEG_FROM_RAD = 180.0 / Math.PI;
    this.RAD_FROM_DEG = Math.PI / 180.0;
    this.HOURS_FROM_RAD = 12.0 / Math.PI;
    this.RAD_FROM_HOURS = Math.PI / 12.0;

    this.CosDeg = function (degrees)
    {
        return Math.cos (this.RAD_FROM_DEG * degrees);
    }

    this.SinDeg = function (degrees)
    {
        return Math.sin (this.RAD_FROM_DEG * degrees);
    }
    
    this.TanDeg = function (degrees) 
    {
        return Math.tan (this.RAD_FROM_DEG * degrees);
    }
    
    this.CosHour = function (hours)
    {
        return Math.cos (this.RAD_FROM_HOURS * hours);
    }

    this.SinHour = function (hours)
    {
        return Math.sin (this.RAD_FROM_HOURS * hours);
    }
    
    this.AtanDeg2 = function (y, x)
    {
        return this.DEG_FROM_RAD * Math.atan2 (y, x);
    }
    
    this.FixHours = function (hours)
    {
        return this.FixCycle (hours, 24.0);
    }
    
    this.FixDegrees = function (degrees)
    {
        return this.FixCycle (degrees, 360.0);
    }
    
    this.FixCycle = function (angle, cycle)
    {
        var fraction = angle / cycle;
        return cycle * (fraction - Math.floor (fraction));
    }
    
    this.Polar = function (x, y, z)
    {
        var rho = (x * x) + (y * y);
        var radius = Math.sqrt (rho + (z * z));
        var phi = Angle.AtanDeg2 (y, x);
        if (phi < 0) {
            phi += 360.0;
        }
        var rho = Math.sqrt (rho);
        var theta = Angle.AtanDeg2 (z, rho);
        return new SphericalCoordinates (phi, theta, radius);
    }
    
    this.DMS = function (x)
    {
        var a = {};
        
        a.negative = (x < 0);
        if (a.negative) {
            x = -x;
        }
        
        a.degrees = Math.floor (x);
        x = 60.0 * (x - a.degrees);
        a.minutes = Math.floor (x);
        x = 60.0 * (x - a.minutes);
        a.seconds = Math.round (10.0 * x) / 10.0;   // Round to the nearest tenth of an arcsecond.
        
        if (a.seconds == 60) {
            a.seconds = 0;
            if (++a.minutes == 60) {
                a.minutes = 0;
                ++a.degrees;
            }
        }
        
        return a;
    }
    
    this.DMM = function (x)
    {
        var a = {};
        a.negative = (x < 0);
        if (a.negative) {
            x = -x;
        }
        
        a.degrees = Math.floor (x);
        x = 60.0 * (x - a.degrees);
        a.minutes = Math.round (100.0 * x) / 100.0;     // Round to nearest hundredth of an arcminute.
        a.seconds = 0.0;        // Fill in just for consistency with Angle.DMS
        
        if (a.minutes >= 60.0) {
            a.minutes -= 60.0;
            ++a.degrees;
        }
        
        return a;
    }
    
    this.SafeArcSinInDegrees = function (z)
    {
        var abs = Math.abs (z);
        if (abs > 1.0) {
            // Guard against blowing up due to slight roundoff errors in Math.Asin ...
            if (abs > 1.00000001) {
                throw "Invalid argument to SafeArcSinInDegrees";
            } else if (z < -1.0) {
                return -90.0;
            } else {
                return +90.0;
            }
        } else {
            return this.DEG_FROM_RAD * Math.asin(z);
        }
    }
}

var Angle = new AngleClass();


//----------------------------------------------------------------------------------------------

function DefaultGeocentricCoordinates (day)     // this function is used for most celestial bodies (but not the Moon)
{
    var bc = this.EclipticCartesianCoordinates (day);               // get this body's heliocentric coordinates
    var ec = Astronomy.Earth.EclipticCartesianCoordinates (day);    // get Earth's heliocentric coordinates
    return bc.Subtract (ec);        // subtract vectors to get the vector from Earth to this body
}

function DefaultEclipticAngularCoordinates (day)
{
    var hc = this.EclipticCartesianCoordinates (day);
    var ac = Angle.Polar (hc.x, hc.y, hc.z);
    return ac;
}

function DefaultHorizontalCoordinates (day, location)
{
    var sky = this.EquatorialCoordinates (day, location);
    return new HorizontalCoordinates (sky, location, day);
}

function Log10 (x)
{
    return Math.log(x) / Math.LN10;
}

//----------------------------------------------------------------------------------------------

function SunClass()
{
    this.Name = "Sun";
    
    this.EclipticCartesianCoordinates = function(day)
    {
        return new CartesianCoordinates (0.0, 0.0, 0.0);
    }

    this.EclipticAngularCoordinates = function(day)
    {
        return new SphericalCoordinates (0.0, 0.0, 0.0);
    }
    
    this.GeocentricCoordinates = DefaultGeocentricCoordinates;
    this.EquatorialCoordinates = EqCoords;
    
    this.DistanceFromEarth = function (day)
    {
        return Astronomy.Earth.EclipticCartesianCoordinates(day).Distance();
    }
    
    this.DistanceFromSun = function (day)
    {
        return 0.0;
    }
    
    this.HorizontalCoordinates = DefaultHorizontalCoordinates;
    
    this.RadiusInMeters = 6.96e+8;
    
    this.VisualMagnitude = function (day)
    {
        var e = this.DistanceFromEarth (day);
        return -26.73 + (5.0 * Log10(e));
    }
}

//----------------------------------------------------------------------------------------------

function EarthClass()
{
    this.Name = "Earth";

    this.EclipticCartesianCoordinates = function(day)
    {
        // We use formulas for finding the Sun as seen from Earth, 
        // then negate the (x,y,z) coordinates obtained to get the Earth's position 
        // from the Sun's perspective.

        // http://www.astro.uio.no/~bgranslo/aares/calculate.html       <== Note error in formula for DS, using sin(RS) where it should say sin(LS)
        // http://www.meteorobs.org/maillist/msg09197.html              <== Correct formulas, more accurate (complex)

        // These formulas use 'd' based on days since 1/Jan/2000 12:00 UTC ("J2000.0"), instead of 0/Jan/2000 0:00 UTC ("day value").
        // Correct by subtracting 1.5 days...
        var d = day - 1.5;
        var T = d / 36525.0;                     // Julian centuries since J2000.0
        var L0 = 280.46645 + (36000.76983 * T) + (0.0003032 * T * T);      // Sun's mean longitude, in degrees
        var M0 = 357.52910 + (35999.05030 * T) - (0.0001559 * T * T) - (0.00000048 * T * T * T);      // Sun's mean anomaly, in degrees

        var C =       // Sun's equation of center in degrees
            (1.914600 - 0.004817 * T - 0.000014 * T * T) * Angle.SinDeg (M0) +
            (0.01993 - 0.000101 * T) * Angle.SinDeg (2 * M0) +
            0.000290 * Angle.SinDeg (3 * M0)
        ;

        var LS = L0 + C;         // true ecliptical longitude of Sun

        var e = 0.016708617 - T * (0.000042037 + (0.0000001236 * T));    // The eccentricity of the Earth's orbit.
        var distanceInAU = (1.000001018 * (1 - e * e)) / (1 + e * Angle.CosDeg(M0 + C));     // distance from Sun to Earth in astronomical units (AU)
        var x = -distanceInAU * Angle.CosDeg (LS);
        var y = -distanceInAU * Angle.SinDeg (LS);
        return new CartesianCoordinates (x, y, 0.0);      // the Earth's center is always on the plane of the ecliptic (z=0), by definition!
    }
    
    this.EclipticAngularCoordinates = DefaultEclipticAngularCoordinates;
    
    this.GeocentricCoordinates = function(day)
    {
        return new CartesianCoordinates (0.0, 0.0, 0.0);
    }
    
    this.EquatorialCoordinates = EqCoords;

    this.DistanceFromEarth = function(day)
    {
        return 0.0;     // included for completeness and consistency of planet interface
    }
    
    this.DistanceFromSun = function (day)
    {
        return this.EclipticCartesianCoordinates(day).Distance();
    }
    
    this.VisualMagnitude = function (day)
    {
        throw "Cannot calculate visual magnitude for Earth!";
    }
}

//----------------------------------------------------------------------------------------------
//  We treat Pluto as a special case.  See http://www.stjarnhimlen.se/comp/ppcomp.html#14

function PlutoClass()
{
    this.Name = "Pluto";
    this.EclipticCartesianCoordinates = function(day)
    {
        var S  =   50.03  +  (0.033459652 * day);
        var P  =  238.95  +  (0.003968789 * day);

        var lonecl = 238.9508  +  (0.00400703 * day) -
            19.799 * Angle.SinDeg(  P)  + 19.848 * Angle.CosDeg(  P) +
             0.897 * Angle.SinDeg(2*P)  -  4.956 * Angle.CosDeg(2*P) +
             0.610 * Angle.SinDeg(3*P)  +  1.211 * Angle.CosDeg(3*P) -
             0.341 * Angle.SinDeg(4*P)  -  0.190 * Angle.CosDeg(4*P) +
             0.128 * Angle.SinDeg(5*P)  -  0.034 * Angle.CosDeg(5*P) -
             0.038 * Angle.SinDeg(6*P)  +  0.031 * Angle.CosDeg(6*P) +
             0.020 * Angle.SinDeg(S-P)  -  0.010 * Angle.CosDeg(S-P)
        ;

        var latecl =  -3.9082 -
            5.453 * Angle.SinDeg(  P)   - 14.975 * Angle.CosDeg(  P) +
            3.527 * Angle.SinDeg(2*P)   +  1.673 * Angle.CosDeg(2*P) -
            1.051 * Angle.SinDeg(3*P)   +  0.328 * Angle.CosDeg(3*P) +
            0.179 * Angle.SinDeg(4*P)   -  0.292 * Angle.CosDeg(4*P) +
            0.019 * Angle.SinDeg(5*P)   +  0.100 * Angle.CosDeg(5*P) -
            0.031 * Angle.SinDeg(6*P)   -  0.026 * Angle.CosDeg(6*P) +
            0.011 * Angle.CosDeg(S-P)
        ;

        var r =  40.72 +
            6.68 * Angle.SinDeg(  P)   + 6.90 * Angle.CosDeg(  P) -
            1.18 * Angle.SinDeg(2*P)   - 0.03 * Angle.CosDeg(2*P) +
            0.15 * Angle.SinDeg(3*P)   - 0.14 * Angle.CosDeg(3*P)
        ;

        var coslon = Angle.CosDeg (lonecl);
        var sinlon = Angle.SinDeg (lonecl);
        var coslat = Angle.CosDeg (latecl);
        var sinlat = Angle.SinDeg (latecl);

        var xp = r * coslon * coslat;
        var yp = r * sinlon * coslat;
        var zp = r * sinlat;

        return new CartesianCoordinates (xp, yp, zp);      // the Earth's center is always on the plane of the ecliptic (z=0), by definition!
    }
    
    this.HorizontalCoordinates = DefaultHorizontalCoordinates;
    this.EclipticAngularCoordinates = DefaultEclipticAngularCoordinates;
    this.GeocentricCoordinates = DefaultGeocentricCoordinates;
    this.EquatorialCoordinates = EqCoords;
    this.DistanceFromEarth = function (day)
    {
        return this.GeocentricCoordinates(day).Distance();
    }

    this.DistanceFromSun = function (day)
    {
        return this.EclipticCartesianCoordinates(day).Distance();
    }
        
    this.VisualMagnitude = function (day)
    {
        var s = this.DistanceFromSun (day);
        var e = this.DistanceFromEarth (day);
        return 14.0 + (5.0 * Log10 ((e * s) / (31.97177 * 31.982)));   // a hack that ignores phase angle, based on data from http://www.heavens-above.com
    }
}

//----------------------------------------------------------------------------------------------

function PlanetPS (     // See  http://www.stjarnhimlen.se/comp/ppcomp.html#4
    name,       // name of the object, e.g. "Mars"
    N0, Nc,     // N0 = longitude of the ascending node (deg).  Nc = rate of change in deg/day
    i0, ic,     // inclination to the ecliptic (deg)
    w0, wc,     // argument of perihelion (deg)
    a0, ac,     // semi-major axis, or mean distance from Sun (AU)
    e0, ec,     // eccentricity (0=circle, 0..1=ellipse, 1=parabola)
    M0, Mc,     // M0 = mean anomaly  (deg) (0 at perihelion; increases uniformly with time).  Mc ("mean motion") = rate of change in deg/day = 360/period
    magBase,
    magPhaseFactor,
    magNonlinearFactor,
    magNonlinearExponent )
{
    this.Name = name;
    this.BodyType = 'planet';
    this.N0 = N0;
    this.Nc = Nc;
    this.i0 = i0;
    this.ic = ic;
    this.w0 = w0;
    this.wc = wc;
    this.a0 = a0;
    this.ac = ac;
    this.e0 = e0;
    this.ec = ec;
    this.M0 = M0;
    this.Mc = Mc;
    this.magBase = magBase;
    this.magPhaseFactor = magPhaseFactor;
    this.magNonlinearFactor = magNonlinearFactor;
    this.magNonlinearExponent = magNonlinearExponent;
}

PlanetPS.prototype.HorizontalCoordinates = DefaultHorizontalCoordinates;
PlanetPS.prototype.EclipticAngularCoordinates = DefaultEclipticAngularCoordinates;

PlanetPS.prototype.MeanAnomaly = function (day)     // day = number of days elapsed since December 31, 1999 0:00 UTC, i.e. JD 2451543.5
{
    return this.M0 + (day * this.Mc);
}

PlanetPS.prototype.NodeLongitude = function (day)
{
    return this.N0 + (day * this.Nc);
}

PlanetPS.prototype.Perihelion = function (day)
{
    return this.w0 + (day * this.wc);
}

PlanetPS.prototype.GeocentricCoordinates = DefaultGeocentricCoordinates;
PlanetPS.prototype.EquatorialCoordinates = EqCoords;

PlanetPS.prototype.EclipticCartesianCoordinates = function (day)
{
    var a = this.a0 + (day * this.ac);
    var e = this.e0 + (day * this.ec);
    var M = this.M0 + (day * this.Mc);
    var N = this.N0 + (day * this.Nc);
    var w = this.w0 + (day * this.wc);
    var i = this.i0 + (day * this.ic);
    var E = EccentricAnomaly (e, M);
    
    // Calculate the body's position in its own orbital plane, and its distance from the thing it is orbiting.
    var xv = a * (Angle.CosDeg(E) - e);
    var yv = a * (Math.sqrt(1.0 - e*e) * Angle.SinDeg(E));
    
    var v = Angle.AtanDeg2 (yv, xv);        // True anomaly in degrees: the angle from perihelion of the body as seen by the Sun.
    var r = Math.sqrt (xv*xv + yv*yv);      // Distance from the Sun to the planet in AU

    var cosN  = Angle.CosDeg (N);
    var sinN  = Angle.SinDeg (N);
    var cosi  = Angle.CosDeg (i);
    var sini  = Angle.SinDeg (i);
    var cosVW = Angle.CosDeg (v + w);
    var sinVW = Angle.SinDeg (v + w);

    // Now we are ready to calculate (unperturbed) ecliptic cartesian heliocentric coordinates.
    var xh = r * (cosN*cosVW - sinN*sinVW*cosi);
    var yh = r * (sinN*cosVW + cosN*sinVW*cosi);
    var zh = r * sinVW * sini;

    return this.Perturb (xh, yh, zh, day);
}

PlanetPS.prototype.Perturb = function (xh, yh, zh, day)
{
    // By default we apply no perturbations.
    // Some planets will override this method so that they can correct for
    // the gravitational influences of other planets.
    return new CartesianCoordinates (xh, yh, zh);
}

PlanetPS.prototype.DistanceFromEarth = function (day)   // returns the distance of this planet from Earth in astronomical units (AU)
{
    return this.GeocentricCoordinates(day).Distance();
}


PlanetPS.prototype.DistanceFromSun = function (day)
{
    return this.EclipticCartesianCoordinates(day).Distance();
}


PlanetPS.prototype.VisualMagnitude = function (day)
{
    var distEarth = this.DistanceFromEarth (day);
    var distSun = this.DistanceFromSun (day);
    var phase = Astronomy.SunEarthPhaseAngle (this, day);
    var mag = this.magBase + (5 * Log10 (distSun * distEarth)) + (this.magPhaseFactor * phase);
    if (this.magNonlinearExponent > 0) {
        mag += this.magNonlinearFactor * Math.pow (phase, this.magNonlinearExponent);
    }
    return mag;
}
    

function EccentricAnomaly (e, M)
{
    var E = M + (e * Angle.SinDeg(M) * (1.0 + (e * Angle.CosDeg(M))));
    for(;;) {
        var F = E - (E - (Angle.DEG_FROM_RAD * e * Angle.SinDeg (E)) - M) / (1 - e * Angle.CosDeg (E));
        var error = Math.abs (F - E);
        E = F;
        if (error < 1.0e-8) {
            break;  // the angle is good enough now for our purposes
        }
    }

    return E;
}


function PerturbMajorPlanet (xh, yh, zh, d)
{
    var ecliptic = Astronomy.EclipticLatLon (xh, yh, zh);
    var lonecl = ecliptic.longitude;
    var latecl = ecliptic.latitude;

    var r = Math.sqrt (xh*xh + yh*yh + zh*zh);

    lonecl += this.PerturbEclipticLongitude (d);
    latecl += this.PerturbEclipticLatitude (d);

    var coslon = Angle.CosDeg (lonecl);
    var sinlon = Angle.SinDeg (lonecl);
    var coslat = Angle.CosDeg (latecl);
    var sinlat = Angle.SinDeg (latecl);

    var xp = r * coslon * coslat;
    var yp = r * sinlon * coslat;
    var zp = r * sinlat;

    return new CartesianCoordinates (xp, yp, zp);
}


function PerturbEclipticLongitude_Jupiter (d)
{
    var Mj = Astronomy.Jupiter.MeanAnomaly (d);
    var Ms = Astronomy.Saturn.MeanAnomaly (d);

    var deltaLong =
        -0.332 * Angle.SinDeg(2*Mj - 5*Ms - 67.6)  -
         0.056 * Angle.SinDeg(2*Mj - 2*Ms + 21  )  +
         0.042 * Angle.SinDeg(3*Mj - 5*Ms + 21  )  -
         0.036 * Angle.SinDeg(  Mj - 2*Ms       )  +
         0.022 * Angle.CosDeg(  Mj -   Ms       )  +
         0.023 * Angle.SinDeg(2*Mj - 3*Ms + 52  )  -
         0.016 * Angle.SinDeg(  Mj - 5*Ms - 69  )
    ;

    return deltaLong;
}


function PerturbEclipticLatitude_Jupiter (d)
{
    return 0.0;
}


function PerturbEclipticLongitude_Saturn (d)
{
    var Mj = Astronomy.Jupiter.MeanAnomaly (d);
    var Ms = Astronomy.Saturn.MeanAnomaly (d);

    var deltaLong =
         0.812 * Angle.SinDeg (2*Mj - 5*Ms - 67.6) -
         0.229 * Angle.CosDeg (2*Mj - 4*Ms -  2.0) +
         0.119 * Angle.SinDeg (  Mj - 2*Ms -  3.0) +
         0.046 * Angle.SinDeg (2*Mj - 6*Ms - 69.0) +
         0.014 * Angle.SinDeg (  Mj - 3*Ms + 32.0)
    ;

    return deltaLong;
}


function PerturbEclipticLatitude_Saturn (d)
{
    var Mj = Astronomy.Jupiter.MeanAnomaly (d);
    var Ms = Astronomy.Saturn.MeanAnomaly (d);

    var deltaLat =
        -0.020 * Angle.CosDeg(2*Mj - 4*Ms - 2)  + 
         0.018 * Angle.SinDeg(2*Mj - 6*Ms - 49)
    ;

    return deltaLat;
}


function PerturbEclipticLongitude_Uranus (d)
{
    var Mj = Astronomy.Jupiter.MeanAnomaly (d);
    var Ms = Astronomy.Saturn.MeanAnomaly (d);
    var Mu = this.MeanAnomaly (d);

    var deltaLong =
        +0.040 * Angle.SinDeg(Ms - 2*Mu + 6)
        +0.035 * Angle.SinDeg(Ms - 3*Mu + 33)
        -0.015 * Angle.SinDeg(Mj - Mu + 20)
    ;

    return deltaLong;
}

function PerturbEclipticLatitude_Uranus (d)
{
    return 0.0;
}

function EqCoords (day, location)
{
    var vectorFromEarth = this.GeocentricCoordinates (day);
    var dx = vectorFromEarth.x;
    var dy = vectorFromEarth.y;
    var dz = vectorFromEarth.z;

    // Now (dx,dy,dz) comprise the vector from the Earth to the object in cartesian ecliptic coordinates.
    // We convert here to equatorial coordinates, using formulas based on the precession of the Earth's axis of rotation.
    var T = Astronomy.CenturiesSinceJ2000 (day);
    var K = 23.4392911 - ((46.8150 * T) - (0.00059 * T * T) + (0.001813 * T * T * T))/3600.0;    // obliquity of ecliptic in degrees.
    var cosK = Angle.CosDeg(K);
    var sinK = Angle.SinDeg(K);

    // Calculate equatorial cartesian coordinates using ecliptic cartesian coordinates...
    var qx = dx;
    var qy = (dy * cosK) - (dz * sinK);
    var qz = (dy * sinK) + (dz * cosK);

    var eq = Angle.Polar (qx, qy, qz);
    eq.longitude /= 15.0;       // convert degrees to sidereal hours
    
    var DEC = eq.latitude;
    var RA  = eq.longitude;     
    var distanceInAU = eq.radius;

    // Determine whether a topocentric correction is warranted.
    // Imagine the angle between the the center of the Earth (geocentric), the remote body,
    // and the observer on the surface of the Earth (topocentric).
    // If this angle is less than 1 arcminute, we won't bother with the extra calculation time.
    // Use approximation to determine whether the parallax matters:
    // asin(RE/r) >= 1"              where RE = radius of Earth, r = distance to CelestialBody.
    // RE/r >= pi / (180 * 3600)
    var parallaxInRadians = 1.0 / (EARTH_RADII_PER_ASTRONOMICAL_UNIT * distanceInAU);
    if (parallaxInRadians >= Math.PI / (180.0 * 3600.0)) {
        // We were using the approximation that parallax = asin(parallax).
        // Now that we know the parallax is significant, go ahead and calculate the exact angle.
        parallaxInRadians = Math.asin (parallaxInRadians);

        // It is easier to calculate topocentric correction in horizontal coordinates...
        var hor = new HorizontalCoordinates (eq, location, day, 0.0);
        
        var gclat;
        var rho;
        
        if (Astronomy.CorrectForOblateEarth) {
            gclat = OblateLatitudeCorrection (location.latitude);
            rho   = OblateRadiusCorrection   (location.latitude);
        } else {
            gclat = location.latitude;
            rho   = 1.0;
        }

        var altitude = hor.altitude - ((Angle.DEG_FROM_RAD * parallaxInRadians) * Angle.CosDeg (hor.altitude));
        var Ls = MeanLongitudeOfSun (day);
        var GMST0 = (Ls + 180.0) / 15.0;
        var utcHours = (day - Math.floor(day)) * 24.0;
        var LST = GMST0 + utcHours + (location.longitude / 15.0);
        var hourAngle = LST - eq.longitude;
        var g = Angle.DEG_FROM_RAD * (Math.atan(Angle.TanDeg(gclat) / Angle.CosHour(hourAngle)));

        var topRA   = eq.longitude  - ((parallaxInRadians * Angle.HOURS_FROM_RAD) * rho * Angle.CosDeg(gclat) * Angle.SinHour(hourAngle) / Angle.CosDeg(eq.latitude));
        var topDEC;
        
        if (Math.abs(g) < 1.0e-6) {
            topDEC = eq.latitude - ((parallaxInRadians * Angle.DEG_FROM_RAD) * rho * Angle.SinDeg(-eq.latitude) * Angle.CosHour(hourAngle));
        } else {
            topDEC = eq.latitude -  (parallaxInRadians * Angle.DEG_FROM_RAD) * rho * Angle.SinDeg(gclat) * Angle.SinDeg(g - eq.latitude) / Angle.SinDeg(g);
        }

        eq = new SphericalCoordinates (topRA, topDEC, distanceInAU);
    }

    return eq;
}

function OblateLatitudeCorrection (latitude)
{
    // This function corrects for the flattening of the Earth due to its rotation.
    // Given a geographic latitude (which is defined based on the tilt of one's local horizon),
    // this function returns geocentric latitude (based on the angle between the equatorial plane and the location).
    // The correction is zero at the equator and at both poles, and is maximized at |latitude| = 45 degrees.
    // See:
    //      http://en.wikipedia.org/wiki/Latitude#Common_.22latitude.22
    //      http://en.wikipedia.org/wiki/Latitude#Geocentric_latitude
    return latitude - (0.1924 * Angle.SinDeg(2.0 * latitude));
}

function OblateRadiusCorrection (latitude)
{
    // Returns the fraction of equatorial Earth radius for sea level at the given geographic latitude.
    // This is due to flattening caused by Earth's rotation.
    // When latitude==0 (i.e. a point on the equator), the value returned is 1.0.
    // The value is minimized at latitude==+90 (North Pole) or latitude==-90 (South Pole).
    return 0.99833 + (0.00167 * Angle.CosDeg(2.0 * latitude));
}

function HorizontalCoordinates (sky, location, day, horizonCorrectionInArcMinutes)
{
    // http://en.wikipedia.org/wiki/Horizontal_coordinate_system
    
    if (horizonCorrectionInArcMinutes == null) {
        horizonCorrectionInArcMinutes = -34.0;
    }

    var GST = GreenwichSiderealTimeInHours (day);
    var LST = GST + (location.longitude / 15.0);
    var hourAngle = Angle.FixHours (LST - sky.longitude);
    
    var sinLat = Angle.SinDeg (location.latitude);
    var cosLat = Angle.CosDeg (location.latitude);

    var sinHourAngle = Angle.SinHour (hourAngle);
    var cosHourAngle = Angle.CosHour (hourAngle);

    var sinDec = Angle.SinDeg (sky.latitude);
    var cosDec = Angle.CosDeg (sky.latitude);

    var altitudeRatio = (sinLat * sinDec) + (cosLat * cosDec * cosHourAngle);
    // Correct for values that are out of range for inverse sine function...
    var absAltitudeRatio = Math.abs (altitudeRatio);
    if (absAltitudeRatio > 1.0) {
        if (absAltitudeRatio > 1.000001) {
            // This is an excessive amount of error: something must be wrong with the formula!
            throw ("Internal error: altitude would be a complex number!");
        } else {
            // Just correct for apparent roundoff without going bezerk.
            this.altitude = (altitudeRatio < 0) ? -90.0 : +90.0;
            this.azimuth = 180.0;    // doesn't really matter what angle we assign: use this value to assist culmination algorithm.
        }
    } else {
        this.altitude = Angle.DEG_FROM_RAD * Math.asin (altitudeRatio);
        var absAltitude = Math.abs (this.altitude);
        var ANGLE_CORRECTION_WINDOW = 6.0;      // I chose 6 degrees as the refraction cutoff, because that is used for Civil Twilight, which should not be corrected.
        if ((absAltitude < ANGLE_CORRECTION_WINDOW) && (horizonCorrectionInArcMinutes != 0.0)) {
            // http://www.jgiesen.de/refract/index.html  For more technical details.
            // Correct for atmospheric lensing making the object appear higher than it would if the Earth had no atmosphere.
            // We want the correction to maximize its effect at the horizon, and vanish to zero as |altitude| approaches ANGLE_CORRECTION_WINDOW.
            // For simplicity we apply a linear interpolation within this range of altitudes.
            var linearFactor = (ANGLE_CORRECTION_WINDOW - absAltitude) / ANGLE_CORRECTION_WINDOW;
            this.altitude -= (horizonCorrectionInArcMinutes / 60.0) * linearFactor;
        }

        this.azimuth = Angle.FixDegrees (Angle.AtanDeg2 (-cosDec * sinHourAngle, (cosLat * sinDec) - (sinLat * cosDec * cosHourAngle)));
    }
}


function GreenwichSiderealTimeInHours (day)
{
    var midnight = Math.floor (day);    // discard fractional part to get same calendar day at midnight UTC
    var T0 = Astronomy.CenturiesSinceJ2000 (midnight);
    var tUT = (day - midnight) * 24.0;     // Greenwich time of day, in hours
    var SG = (6.6974 + 2400.0513 * T0) + (366.2422 / 365.2422) * tUT;
    SG = Angle.FixHours (SG);
    return SG;
}


//----------------------------------------------------------------------------------------------

function CreateAsteroid (
    name,       // name of asteroid
    epochJD,    // epoch of orbital elements as Julian Date
    Nx,         // longitude of the ascending node at epoch, in degrees
    i,          // inclination to the ecliptic, in degrees
    w,          // argument of perihelion, in degrees
    a,          // semi-major axis in AU
    e,          // eccentricity
    Mx,         // mean anomaly at epoch, in degrees
    T,          // orbital period, in days
    amag )      // absolute magnitude
{
    var day = epochJD - 2451543.5;              // convert Julian Date to "0.0 January 2000" standard epoch day value.
    var Mc = 360.0 / T;                         // "mean motion": how many degrees per day the body orbits around the Sun, on average.
    var M0 = Angle.FixDegrees (Mx - Mc*day);    // work backwards to figure out mean anomoly at standard epoch.

    var N0 = Nx;    //!!! FIXFIXFIX
    var Nc = 0.0;   //!!! FIXFIXFIX

    var p = new PlanetPS (
        name,
        N0, Nc,
        i, 0.0,
        w, 0.0,
        a, 0.0,
        e, 0.0,
        M0, Mc,
        amag,
        0.0,
        0.0,
        0.0
    );
    
    p.BodyType = 'asteroid';        // change body type from 'planet' to 'asteroid'
    
    return p;
}

function CreateComet(/*see arguments for CreateAsteroid*/) {
    // Pass along the same arguments to "CreateAsteroid".
    var comet = CreateAsteroid.apply(this, arguments);
    comet.BodyType = 'comet';       // change body type from 'asteroid' to 'comet'.
    return comet;
}

function CreateMinor(/*see arguments for CreateAsteroid*/) {
    // Pass along the same arguments to "CreateAsteroid".
    var comet = CreateAsteroid.apply(this, arguments);
    comet.BodyType = 'minor';       // indicate that this is a minor asteroid
    return comet;
}

function CreatePlanetJPL(       //  Designed for use with JPL orbital data
    name,       // name of the object, e.g. "Mars"
    ja0, jac,   // semi-major axis (AU), rate (AU/century)
    je0, jec,   // eccentricity (1), (1/century)
    jI0, jIc,   // inclination (deg), (deg/century)
    jL0, jLc,   // mean longitude (deg), (deg/century)
    ju0, juc,   // longitude of perihelion (deg), (deg/century)
    jQ0, jQc,   // longitude of ascending node (deg), (deg/century)
    magBase,
    magPhaseFactor,
    magNonlinearFactor,
    magNonlinearExponent )
{
    // Convert units and epoch from JPL format to PlanetPS format.
    // http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
    // http://ssd.jpl.nasa.gov/txt/p_elem_t1.txt
    var dpc = 36525.0;      // days per century
    var cofs = 1.5/dpc;     // centuries that JPL epoch (J2000 = 1/1/2000 12UT) is ahead of PS (12/31/1999 0UT)
    
    var N0 = jQ0 - cofs*jQc;
    var Nc = jQc/dpc;
    
    var i0 = jI0 - cofs*jIc;
    var ic = jIc/dpc;
    
    // w = ju - jQ
    var w0 = Angle.FixDegrees((ju0 - cofs*juc) - (jQ0 - cofs*jQc));
    var wc = (juc - jQc)/dpc;
    
    var a0 = ja0 - cofs*jac;
    var ac = jac/dpc;
    
    var e0 = je0 - cofs*jec;
    var ec = jec/dpc;
    
    // M = jL - ju
    var M0 = Angle.FixDegrees((jL0 - cofs*jLc) - (ju0 - cofs*juc));
    var Mc = (jLc - juc)/dpc;

    return new PlanetPS(
        name,       // name of the object, e.g. "Mars"
        N0, Nc,     // N0 = longitude of the ascending node (deg).  Nc = rate of change in deg/day
        i0, ic,     // inclination to the ecliptic (deg)
        w0, wc,     // argument of perihelion (deg)
        a0, ac,     // semi-major axis, or mean distance from Sun (AU)
        e0, ec,     // eccentricity (0=circle, 0..1=ellipse, 1=parabola)
        M0, Mc,     // M0 = mean anomaly  (deg) (0 at perihelion; increases uniformly with time).  Mc ("mean motion") = rate of change in deg/day = 360/period
        magBase,
        magPhaseFactor,
        magNonlinearFactor,
        magNonlinearExponent
    );
}

function CreateJupiter() {
    var planet = CreatePlanetJPL(
        "Jupiter",
        5.20288700,  -0.00011607,    // ja0, jac = semi-major axis (AU), rate (AU/century)
        0.04838624,  -0.00013253,    // je0, jec = eccentricity (1), (1/century)
        1.30439695,   0.00183714,    // jI0, jIc = inclination (deg), (deg/century)
        34.39644051, 3034.74612775,  // jL0, jLc = mean longitude (deg), (deg/century)
        14.72847983,  0.21252668,    // ju0, juc = longitude of perihelion (deg), (deg/century)
        100.47390909, 0.20469106,    // jQ0, jQc = longitude of ascending node (deg), (deg/century)
        -9.25, 0.014, 0, 0      // visual magnitude elements)       
    );
    
    planet.Perturb = PerturbMajorPlanet;   // override the do-nothing method from PlanetPS
    planet.PerturbEclipticLongitude = PerturbEclipticLongitude_Jupiter;     // method called by PerturbMajorPlanet
    planet.PerturbEclipticLatitude  = PerturbEclipticLatitude_Jupiter;      // method called by PerturbMajorPlanet
    
    return planet;
}

function CreateSaturnJPL()
{
    var planet = CreatePlanetJPL (
        "SaturnJPL",
         9.53667594,   -0.00125060,    // ja0, jac = semi-major axis (AU), rate (AU/century)
         0.05386179,   -0.00050991,    // je0, jec = eccentricity (1), (1/century)
         2.48599187,    0.00193609,    // jI0, jIc = inclination (deg), (deg/century)
        49.95424423, 1222.49362201,    // jL0, jLc = mean longitude (deg), (deg/century)
        92.59887831,   -0.41897216,    // ju0, juc = longitude of perihelion (deg), (deg/century)
       113.66242448,   -0.28867794,    // jQ0, jQc = longitude of ascending node (deg), (deg/century)
        -9.0, 0.044, 0, 0     // visual magnitude elements
    );
    
    planet.Perturb = PerturbMajorPlanet;   // override the do-nothing method from PlanetPS
    planet.PerturbEclipticLongitude = PerturbEclipticLongitude_Saturn;     // method called by PerturbMajorPlanet
    planet.PerturbEclipticLatitude  = PerturbEclipticLatitude_Saturn;      // method called by PerturbMajorPlanet
    planet.__BaseVisualMagnitude = planet.VisualMagnitude;    // we override the base function to also include a factor for the rings...
    
    planet.VisualMagnitude = function (day)
    {
        var planetMagnitude = this.__BaseVisualMagnitude(day);  // get magnitude of the planet body itself.
        
        var Ir = 28.06;    // inclination of Saturn's rings to ecliptic, in degrees
        var cosIr = Angle.CosDeg (Ir);
        var sinIr = Angle.SinDeg (Ir);

        var gc = this.GeocentricCoordinates (day);
        var Los = Angle.FixDegrees (Angle.AtanDeg2 (gc.y, gc.x));
        var Las = Angle.FixDegrees (Angle.AtanDeg2 (gc.z, Math.sqrt(gc.x*gc.x + gc.y*gc.y)));
        var sinLas = Angle.SinDeg (Las);
        var cosLas = Angle.CosDeg (Las);
        var Nr = 169.51 + (3.82e-5 * day);     // ascending node of the plane of Saturn's rings
        var sinLosMinusNr = Angle.SinDeg (Los - Nr);

        var B = Math.asin (sinLas*cosIr - cosLas*sinIr*sinLosMinusNr);
        var sinB = Math.abs (Math.sin (B));     // ??? can we get rid of doing both Asin and Sin?

        var ringMagnitude = (-2.6 * Math.abs (sinB)) + (1.2 * sinB * sinB);
        return planetMagnitude + ringMagnitude;
    }
    
    return planet;
}


function CreateSaturn()
{
    var planet = new PlanetPS (
        "Saturn", 
        113.6634, 2.3898e-5, 2.4886, -1.081e-7, 339.3939, 2.97661e-5, 9.55475, 0.0, 0.055546, -9.499e-9, 316.9670, 0.0334442282,    // orbital elements
        -9.0, 0.044, 0, 0     // visual magnitude elements
    );
    
    planet.Perturb = PerturbMajorPlanet;   // override the do-nothing method from PlanetPS
    planet.PerturbEclipticLongitude = PerturbEclipticLongitude_Saturn;     // method called by PerturbMajorPlanet
    planet.PerturbEclipticLatitude  = PerturbEclipticLatitude_Saturn;      // method called by PerturbMajorPlanet
    planet.__BaseVisualMagnitude = planet.VisualMagnitude;    // we override the base function to also include a factor for the rings...
    
    planet.VisualMagnitude = function (day)
    {
        var planetMagnitude = this.__BaseVisualMagnitude(day);  // get magnitude of the planet body itself.
        
        var Ir = 28.06;    // inclination of Saturn's rings to ecliptic, in degrees
        var cosIr = Angle.CosDeg (Ir);
        var sinIr = Angle.SinDeg (Ir);

        var gc = this.GeocentricCoordinates (day);
        var Los = Angle.FixDegrees (Angle.AtanDeg2 (gc.y, gc.x));
        var Las = Angle.FixDegrees (Angle.AtanDeg2 (gc.z, Math.sqrt(gc.x*gc.x + gc.y*gc.y)));
        var sinLas = Angle.SinDeg (Las);
        var cosLas = Angle.CosDeg (Las);
        var Nr = 169.51 + (3.82e-5 * day);     // ascending node of the plane of Saturn's rings
        var sinLosMinusNr = Angle.SinDeg (Los - Nr);

        var B = Math.asin (sinLas*cosIr - cosLas*sinIr*sinLosMinusNr);
        var sinB = Math.abs (Math.sin (B));     // ??? can we get rid of doing both Asin and Sin?

        var ringMagnitude = (-2.6 * Math.abs (sinB)) + (1.2 * sinB * sinB);
        return planetMagnitude + ringMagnitude;
    }
    
    return planet;
}

function CreateUranus()
{
    var planet = CreatePlanetJPL(
        "Uranus", 
        19.18916464,   -0.00196176,     // ja0, jac = semi-major axis (AU), rate (AU/century)
         0.04725744,   -0.00004397,     // je0, jec = eccentricity (1), (1/century)
         0.77263783,   -0.00242939,     // jI0, jIc = inclination (deg), (deg/century)
       313.23810451,  428.48202785,     // jL0, jLc = mean longitude (deg), (deg/century)
       170.95427630,    0.40805281,     // ju0, juc = longitude of perihelion (deg), (deg/century)
        74.01692503,    0.04240589,     // jQ0, jQc = longitude of ascending node (deg), (deg/century)
        -7.15, 0.001, 0, 0      // visual magnitude elements
    );
    
    planet.Perturb = PerturbMajorPlanet;   // override the do-nothing method from PlanetPS
    planet.PerturbEclipticLongitude = PerturbEclipticLongitude_Uranus;     // method called by PerturbMajorPlanet
    planet.PerturbEclipticLatitude  = PerturbEclipticLatitude_Uranus;      // method called by PerturbMajorPlanet
    
    return planet;
}

//----------------------------------------------------------------------------------------------

function MeanAnomalyOfSun (d)
{
    return 356.0470 + (0.9856002585 * d);
}

function SunArgumentOfPerihelion (d)
{
    return 282.9404 + (4.70935e-5 * d);
}

function MeanLongitudeOfSun (d)
{
    return MeanAnomalyOfSun(d) + SunArgumentOfPerihelion(d);
}

function CreateMoon()
{
    var moon = new PlanetPS (
        "Moon",
        125.1228, -0.0529538083, 5.1454, 0, 318.0634, 0.1643573223, 60.2666 / EARTH_RADII_PER_ASTRONOMICAL_UNIT, 0, 0.054900, 0, 115.3654, 13.0649929509,
        0.23, 0.026, 4.0e-9, 4
    );
    
    moon.EclipticCartesianCoordinates = function(day)
    {
        var mc = this.GeocentricCoordinates (day);
        var ec = Astronomy.Earth.EclipticCartesianCoordinates (day);
        return ec.Add (mc);
    }
    
    moon.GeocentricCoordinates = PlanetPS.prototype.EclipticCartesianCoordinates;   // Moon uses same formulas as planet, but results in geocentric values!
    
    moon.Perturb = function (xh, yh, zh, d)
    {
        var Ms = MeanAnomalyOfSun (d);               // mean anomaly of Sun
        var ws = SunArgumentOfPerihelion (d);        // Sun's argument of perihelion
        var Ls = Ms + ws;                            // mean longitude of Sun

        var Mm = this.MeanAnomaly (d);               // Moon's mean anomaly
        var Nm = this.NodeLongitude (d);             // longitude of Moon's node
        var wm = this.Perihelion (d);                // Moon's argument of perihelion
        var Lm = Mm + wm + Nm;                       // Mean longitude of the Moon

        var D = Lm - Ls;                             // mean elongation of the Moon
        var F = Lm - Nm;                             // argument of latitude for the Moon

        var deltaLong =
            -1.274 * Angle.SinDeg(Mm - 2*D)      +   // the Evection
             0.658 * Angle.SinDeg(2*D)           -   // the Variation
             0.186 * Angle.SinDeg(Ms)            -   // the Yearly Equation
             0.059 * Angle.SinDeg(2*Mm - 2*D)    -
             0.057 * Angle.SinDeg(Mm - 2*D + Ms) +
             0.053 * Angle.SinDeg(Mm + 2*D)      +
             0.046 * Angle.SinDeg(2*D - Ms)      +
             0.041 * Angle.SinDeg(Mm - Ms)       -
             0.035 * Angle.SinDeg(D)             -   // the Parallactic Equation
             0.031 * Angle.SinDeg(Mm + Ms)       -
             0.015 * Angle.SinDeg(2*F - 2*D)     +
             0.011 * Angle.SinDeg(Mm - 4*D)
        ;

        var deltaLat =
            -0.173 * Angle.SinDeg(F - 2*D)       -
             0.055 * Angle.SinDeg(Mm - F - 2*D)  -
             0.046 * Angle.SinDeg(Mm + F - 2*D)  +
             0.033 * Angle.SinDeg(F + 2*D)       +
             0.017 * Angle.SinDeg(2*Mm + F)
        ;

        var deltaRadius =
            -0.58 * Angle.CosDeg (Mm - 2*D)   -
             0.46 * Angle.CosDeg (2*D)
        ;

        var ecliptic = Astronomy.EclipticLatLon (xh, yh, zh);
        var lonecl = ecliptic.longitude;
        var latecl = ecliptic.latitude;

        var r = Math.sqrt (xh*xh + yh*yh + zh*zh);

        lonecl += deltaLong;
        latecl += deltaLat;
        r += deltaRadius / EARTH_RADII_PER_ASTRONOMICAL_UNIT;

        var coslon = Angle.CosDeg (lonecl);
        var sinlon = Angle.SinDeg (lonecl);
        var coslat = Angle.CosDeg (latecl);
        var sinlat = Angle.SinDeg (latecl);

        var xp = r * coslon * coslat;
        var yp = r * sinlon * coslat;
        var zp = r * sinlat;

        return new CartesianCoordinates (xp, yp, zp);
    }
    
    moon.RadiusInMeters = 1.7374e+6;
    
    return moon;
}

//----------------------------------------------------------------------------------------------


function AstronomyClass()
{
    this.DayValue = function(utc)
    {
        // http://www.elated.com/articles/working-with-dates/
        // Return number of days since 0/Jan/2000 00:00 UTC...
        return 1.0 + (utc.getTime() - Date.UTC(2000, 0, 1)) / (3600.0 * 24.0 * 1000.0);
    }
    
    this.DaysSinceJ2000 = function(day)
    {
        return day - 1.5;
    }
    
    this.CenturiesSinceJ2000 = function(day)
    {
        return this.DaysSinceJ2000(day) / 36525.0;
    }
    
    this.CurrentDayValue = function(utc)
    {
        return this.DayValue(new Date());
    }
    
    this.DayValueToDate = function (day)
    {
        var date = new Date();
        var M = 3600.0 * 24.0 * 1000.0;
        var T0 = Date.UTC (2000, 0, 1);
        date.setTime (M*(day - 1.0) + T0);
        return date;
    }
    
    this.EclipticLatLon = function (x, y, z)
    {
        return { 'longitude': Angle.AtanDeg2(y, x), 'latitude': Angle.AtanDeg2(z, Math.sqrt(x*x + y*y)) };
    }
    
    this.NextRiseTime = function (body, day, location)
    {
        return Astronomy_FindNextTransition (body, day, day + 1.1, location, Astronomy_RiseCondition, 3);
    }
    
    this.NextSetTime = function (body, day, location)
    {
        return Astronomy_FindNextTransition (body, day, day + 1.1, location, Astronomy_SetCondition, 3);
    }
    
    this.NextCulmTime = function (body, day, location)
    {
        return Astronomy_FindNextTransition (body, day, day + 1.1, location, Astronomy_CulminateCondition, 6);
    }

    this.NextNewMoon = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (this.Moon, day, day + numDays, null, Astronomy_RelativeLongitudeCondition, numIntervals, 0.0);
    }
    
    this.NextMoonFirstQuarter = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (this.Moon, day, day + numDays, null, Astronomy_RelativeLongitudeCondition, numIntervals, 90.0);
    }

    this.NextFullMoon = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (this.Moon, day, day + numDays, null, Astronomy_RelativeLongitudeCondition, numIntervals, 180.0);
    }

    this.NextMoonThirdQuarter = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (this.Moon, day, day + numDays, null, Astronomy_RelativeLongitudeCondition, numIntervals, 270.0);
    }
    
    this.NextOpposition = function (body, day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (body, day, day + numDays, null, Astronomy_RelativeLongitudeCondition, numIntervals, 180.0);
    }
    
    this.NextConjunction = function (body, day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (body, day, day + numDays, null, Astronomy_RelativeLongitudeCondition, numIntervals, 0.0);
    }
    
    this.NextMaxSunAngle = function (body, day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (body, day, day + numDays, null, Astronomy_MaxSunAngleCondition, numIntervals);
    }
    
    this.NextVernalEquinox = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (this.Sun, day, day + numDays, null, Astronomy_VernalEquinoxCondition, numIntervals);
    }

    this.NextAutumnalEquinox = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (this.Sun, day, day + numDays, null, Astronomy_AutumnalEquinoxCondition, numIntervals);
    }

    this.NextNorthernSolstice = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (this.Sun, day, day + numDays, null, Astronomy_NorthernSolsticeCondition, numIntervals);
    }

    this.NextSouthernSolstice = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (this.Sun, day, day + numDays, null, Astronomy_SouthernSolsticeCondition, numIntervals);
    }
    
    this.NextPeakVisualMagnitude = function (body, day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (body, day, day + numDays, null, Astronomy_PeakVisualMagnitudeCondition, numIntervals);
    }

    this.NextMinAngleWithOtherBody = function (body, day, numDays, numIntervals, otherBody)
    {
        return Astronomy_FindNextTransition (body, day, day + numDays, null, Astronomy_MinAngleWithOtherBodyCondition, numIntervals, otherBody);
    }

    this.NextMoonPerigee = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (null, day, day + numDays, null, Astronomy_MoonPerigee, numIntervals, null);
    }

    this.NextMoonApogee = function (day, numDays, numIntervals)
    {
        return Astronomy_FindNextTransition (null, day, day + numDays, null, Astronomy_MoonApogee, numIntervals, null);
    }
    
    this.SunEarthPhaseAngle = function (body, day)
    {
        if (body.Name == "Sun") {
            throw "Cannot calculate SunEarthPhaseAngle of the Sun";
        } else {
            var h = body.EclipticCartesianCoordinates (day);        // a.k.a. Heliocentric Coordinates
            var g = body.GeocentricCoordinates (day);
            
            // We take advantage of the fact that when two lines cross, opposite angles are equal.
            // In other words, we want the angle between (-h) and (-g), which is the same as the angle between h and g,
            // because (-h) dot (-g) = h dot g :
            // ((-hx)*(-gx) = hx*gx, etc.

            return AngleBetweenVectorsInDegrees (h, g);
        }
    }
    
    this.AngleWithSunInDegrees = function (body, day)
    {
        if (body.Name == "Sun") {
            return 0.0;
        } else {
            var s = Astronomy.Sun.GeocentricCoordinates (day);
            var b = body.GeocentricCoordinates (day);
            return AngleBetweenVectorsInDegrees (s, b);
        }
    }

    this.AngleBetweenBodiesInDegrees = function (body1, body2, day)
    {
        var a = body1.GeocentricCoordinates (day);
        var b = body2.GeocentricCoordinates (day);
        return AngleBetweenVectorsInDegrees (a, b);
    }
    
    this.RelativeLongitude = function (body, day)
    {
        // Returns the relative longitude between the body and the Sun as seen from Earth's center.
        var bc = body.GeocentricCoordinates (day);
        var sc = this.Sun.GeocentricCoordinates (day);
        
        var blon = Angle.AtanDeg2 (bc.y, bc.x);
        var slon = Angle.AtanDeg2 (sc.y, sc.x);
        var rlon = Angle.FixDegrees (blon - slon);
        return rlon;
    }
    
    this.PrecessionRotationMatrix = function (EPOCH1, EPOCH2)
    {
        var CDR = Math.PI / 180.0;
        var CSR = CDR / 3600.0;
        var T = 0.001 * (EPOCH2 - EPOCH1);
        var ST = 0.001 * (EPOCH1 - 1900.0);
        var A = CSR * T * (23042.53 + ST * (139.75 + 0.06 * ST) + T * (30.23 - 0.27 * ST + 18.0 * T));
        var B = CSR * T * T * (79.27 + 0.66 * ST + 0.32 * T) + A;
        var C = CSR * T * (20046.85 - ST * (85.33 + 0.37 * ST) + T * (-42.67 - 0.37 * ST - 41.8 * T));
        var SINA = Math.sin (A);
        var SINB = Math.sin (B);
        var SINC = Math.sin (C);
        var COSA = Math.cos (A);
        var COSB = Math.cos (B);
        var COSC = Math.cos (C);

        var r = [[], [], []];
        r[0][0] = COSA * COSB * COSC - SINA * SINB;
        r[0][1] = -COSA * SINB - SINA * COSB * COSC;
        r[0][2] = -COSB * SINC;
        r[1][0] = SINA * COSB + COSA * SINB * COSC;
        r[1][1] = COSA * COSB - SINA * SINB * COSC;
        r[1][2] = -SINB * SINC;
        r[2][0] = COSA * SINC;
        r[2][1] = -SINA * SINC;
        r[2][2] = COSC;

        return r;
    }
    
    // We are taking a short-cut for computational efficiency:
    // We calculate the rotation matrix for the current date a single time,
    // then keep using it.  This will be plenty accurate for anything
    // less than a decade before or after the current date.
    //
    this.ConstellationRotationMatrix = this.PrecessionRotationMatrix (new Date().getFullYear(), 1875);
    
    /*
        The PrecessEquatorialCoordinates() function was translated from C to C# by Don Cross, 
        then to Javascript by Don Cross, based on
        code with original comments thus:

        This program is a translation with a few adaptations of the 
        Fortran program.f, made by FO @ CDS  (francois@simbad.u-strasbg.fr)  
        in November 1996.
    */
    this.PrecessEquatorialCoordinates = function (eq, matrix)
    {
        // Don says: in the original C code, the RA1 and DEC1 passed in were always in radians.
        // I want to pass in hours and degrees, respectively, so convert to radians here...
        var RA1  = Angle.RAD_FROM_HOURS * eq.longitude;
        var DEC1 = Angle.RAD_FROM_DEG   * eq.latitude;

        /* Compute input direction cosines */
        var A = Math.cos (DEC1);

        var x1 = [
            A * Math.cos (RA1),
            A * Math.sin (RA1),
            Math.sin (DEC1)
        ];

        /* Perform the rotation to get the direction cosines at epoch2 */
        var x2 = [ 0.0, 0.0, 0.0 ];
        for (var i = 0; i < 3; i++) {
            for (var j = 0; j < 3; j++) {
                x2[i] += matrix[i][j] * x1[j];
            }
        }

        // Don says: for my purposes, I always want RA in hours and DEC in degrees...
        var RA2  = Angle.HOURS_FROM_RAD * Math.atan2 (x2[1], x2[0]);
        if (RA2 < 0) {
            RA2 += 24.0;
        }
        
        var DEC2 = Angle.SafeArcSinInDegrees (x2[2]);

        return new SphericalCoordinates (RA2, DEC2, eq.radius);
    }
    
    this.FindConstellation = function (eq)
    {
        var eqOld = this.PrecessEquatorialCoordinates (eq, this.ConstellationRotationMatrix);
        var ra  = eqOld.longitude;
        var dec = eqOld.latitude;
        for (var i=0; i < ConstellationEdgeArray.length; ++i) {
            var e = ConstellationEdgeArray[i];
            if ((e.DecBottom <= dec) && (e.RaRight > ra) && (e.RaLeft <= ra)) {
                var c = ConstellationByConciseName[e.ConciseName];
                return { 'ConciseName':e.ConciseName, 'FullName':c.FullName, 'PossessiveName':c.PossessiveName };
            }
        }
        return null;    // not found!
    }
    
    this.CorrectForOblateEarth = true;      // Set to false if oblateness formulas turn out to have problems.
    
    this.Sun     = new SunClass();
    this.Moon    = CreateMoon();
    
    this.Mercury = new PlanetPS ("Mercury", 48.3313, 3.24587e-5, 7.0047, 5.0e-8, 29.1241, 1.01444e-5, 0.387098, 0.0, 0.205635, 5.59e-10, 168.6562, 4.0923344368, -0.36, 0.027, 2.2e-13, 6);
    this.Venus   = new PlanetPS ("Venus", 76.6799, 2.46590e-5, 3.3946, 2.75e-8, 54.8910, 1.38374e-5, 0.723330, 0.0, 0.006773, -1.302e-9, 48.0052, 1.6021302244, -4.34, 0.013, 4.2e-7, 3);
    this.Earth   = new EarthClass();
    this.Mars    = new PlanetPS ("Mars", 49.5574, 2.11081e-5, 1.8497, -1.78e-8, 286.5016, 2.92961e-5, 1.523688, 0.0, 0.093405, 2.516e-9, 18.6021, 0.5240207766, -1.51, 0.016, 0, 0);

    // Asteroids and comets    
    //                             name       epochJD     Nx=long asc node   i=inclination       w=arg of peri       a=semi-major axis  e=eccentricity       Mx=mean anom@epoch  T=orbital period   abs mag
    this.Ceres   = CreateAsteroid ("Ceres",   2457000.5,  80.32926547452543, 10.59338616262872,  72.52203270043788,  2.76750591440571,  0.07582276595896797,  95.98917578768719, 1681.633163085685,  4.18);   // http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=1
    this.Pallas  = CreateAsteroid ("Pallas",  2457000.5, 173.0962475866183,  34.84099788068097,  309.9303277309815,  2.771606108489468, 0.2312736282433415,   78.22870368538561, 1685.371678425934,  5.15);   // http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=2
    this.Juno    = CreateAsteroid ("Juno",    2457000.5, 169.8711811946533,  12.98166179450518,  248.4099652667037,  2.670700240558908, 0.2554482627817375,   33.07715332505229, 1594.175554183099,  5.91);   // http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=3
    this.Vesta   = CreateAsteroid ("Vesta",   2457000.5, 103.8513672233257,   7.140427316529066, 151.1985275178533,  2.361793227026224, 0.08874017002173754,  20.86384148999364, 1325.748779938602,  3.67);   // http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=4
    this.Ida     = CreateMinor    ("Ida",     2457000.5, 324.0275757719688,   1.132199172379656, 110.4961991586997,  2.862142759005513, 0.04143441719730059, 276.9785940904641,  1768.623391686735,  9.94);   // http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=243
    this.Gaspra  = CreateMinor    ("Gaspra",  2457000.5, 253.1678768421052,   4.102608826356864, 129.4892377057742,  2.209868867807942, 0.1734277243451101,  293.2411112895537,  1199.908645279267, 11.46);   // http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=951
    this.Comet_9P  =  CreateComet ("9P/T1",   2456717.5,  68.88181348869867, 10.50258430337173,  179.3031328527078,  3.140094842040291, 0.5115944387599365,  203.237608505211,   2032.415859334781, 13.10);   // comet 9P/Tempel 1 : http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=9P  
    this.Comet_19P =  CreateComet ("19P/B",   2456870.5,  75.3828362947892,  30.36856213315691,  353.4432880799574,  3.602573738021955, 0.6255352567887734,  316.6290089712147,  2497.570436475743,  8.90);   // comet 19P/Borrelly : http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=19P
    this.Comet_67P =  CreateComet ("67P/C-G", 2456981.5,  50.14210951437195,  7.040200902346087,  12.78560606538363, 3.462817302992186, 0.6409739314162571,  319.3033467788339,  2353.654952366357, 13.55);   // comet 67P/Churyumov-Gerasimenko : http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=67P
    this.Comet_81P =  CreateComet ("81P/W2",  2456978.5, 136.1250968664513,   3.239043597517584,  41.68675738457102, 3.448547391975779, 0.5380434533956898,  265.9016154674808,  2339.121199169129,  8.60);   // comet 81P/Wild 2 : http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=81P 

    this.Jupiter = CreateJupiter();
    this.Saturn  = CreateSaturn();
    //this.SaturnJPL = CreateSaturnJPL();
    this.Uranus  = CreateUranus();
    this.Neptune = new PlanetPS ("Neptune", 131.7806, 3.0173e-5, 1.7700, -2.55e-7, 272.8461, -6.027e-6, 30.05826, 3.313e-8, 0.008606, 2.15e-9, 260.2471, 0.005995147, -6.90, 0.001, 0, 0);
    this.Pluto   = new PlutoClass();
    
    this.Body = [
        this.Sun, 
        this.Mercury, 
        this.Venus, 
        this.Earth, 
        this.Moon, 
        this.Mars, 
        this.Ceres,
        this.Pallas,
        this.Juno,
        this.Vesta, 
        this.Ida,
        this.Gaspra,
        this.Comet_9P,
        this.Comet_19P,
        this.Comet_67P,
        this.Comet_81P,
        this.Jupiter, 
        this.Saturn, 
        //this.SaturnJPL,       // not much better than existing Saturn... not ready for publish
        this.Uranus, 
        this.Neptune, 
        this.Pluto
    ];
    
    this.IsBodyName = function (name)
    {
        return (this[name] != null) && (this[name].Name == name);
    }
    
    this.PlanetNumber = {
        'Mercury':      1,
        'Venus':        2,
        'Earth':        3,
        'Mars':         4,
        'Jupiter':      5,
        'Saturn':       6,
        'Uranus':       7,
        'Neptune':      8,
        'Pluto':        9
    };
}

function Astronomy_AngularRadius (body, day, location)
{
    if (body.RadiusInMeters == null) {
        return 0.0;
    } else {
        // FIXFIXFIX:  use location to figure out distance from body to observer, not just center of Earth!
        var distance = body.DistanceFromEarth (day) * METERS_PER_ASTRONOMICAL_UNIT;        
        return Angle.DEG_FROM_RAD * (body.RadiusInMeters / distance);
    }
}

function Astronomy_UpperLimbAltitude (body, day, location)
{
    var altitude = body.HorizontalCoordinates(day,location).altitude;
    var radius = Astronomy_AngularRadius (body, day, location);
    return altitude + radius;
}

function Astronomy_RiseCondition (body, day1, day2, location)
{
    // Return true if upper limb of body rises in the range of time t : day1 <= t <= day2.
    var altitude1 = Astronomy_UpperLimbAltitude (body, day1, location);
    var altitude2 = Astronomy_UpperLimbAltitude (body, day2, location);
    return (altitude1 <= 0.0) && (altitude2 >= 0.0);
}

function Astronomy_SetCondition (body, day1, day2, location)
{
    // Return true if upper limb of body sets in the range of time t : day1 <= t <= day2.
    var altitude1 = Astronomy_UpperLimbAltitude (body, day1, location);
    var altitude2 = Astronomy_UpperLimbAltitude (body, day2, location);
    return (altitude1 >= 0.0) && (altitude2 <= 0.0);
}

function Astronomy_CulminateCondition (body, day1, day2, location)
{
    // Return true if center of body crosses the meridian between day1 and day2.
    var culm = false;
    var hc1 = body.HorizontalCoordinates (day1, location);
    // Figure out if the object crosses the meridian between day1 and day2.
    // This happens if the object moves from an eastern azimuth to a western azimuth.
    // An azimuth is eastern if sin(azimuth) > 0, or western if sin(azimuth) < 0.
    // The altitude at both times must be above the horizon (altitude > 0).
    if (hc1.altitude > 0) {
        var sin1 = Angle.SinDeg (hc1.azimuth);
        if (sin1 >= 0) {
            var hc2 = body.HorizontalCoordinates (day2, location);
            if (hc2.altitude > 0) {
                var sin2 = Angle.SinDeg (hc2.azimuth);
                if (sin2 <= 0) {
                    culm = true;
                }
            }
        }
    }
    return culm;
}

function Astronomy_MaxSunAngleCondition (body, day1, day2)
{
    // Calculate slopes using a teensy weensy increment.
    // We have found the maximum if the angle is increasing at day1 and decreasing at day2.
    
    var day1b = day1 + 1.0e-7;
    var day2a = day2 - 1.0e-7;
    
    var sa1  = Astronomy.AngleWithSunInDegrees (body, day1 );
    var sa1b = Astronomy.AngleWithSunInDegrees (body, day1b);
    var sa2a = Astronomy.AngleWithSunInDegrees (body, day2a);
    var sa2  = Astronomy.AngleWithSunInDegrees (body, day2 );
    
    return (sa1 <= sa1b) && (sa2a >= sa2);
}

function Astronomy_MoonApogee (body, day1, day2, location, otherBody)
{
    // The parameters 'body', 'location', 'otherBody' are all unused.
    // Return true if Moon is at greatest distance from Earth.
    
    var day1b = day1 + 1.0e-6;
    var day2a = day2 - 1.0e-6;

    var dist1  = Astronomy.Moon.DistanceFromEarth(day1);
    var dist1b = Astronomy.Moon.DistanceFromEarth(day1b);
    var dist2a = Astronomy.Moon.DistanceFromEarth(day2a);
    var dist2  = Astronomy.Moon.DistanceFromEarth(day2);

    return (dist1b >= dist1) && (dist2a >= dist2);
}

function Astronomy_MoonPerigee (body, day1, day2, location, otherBody)
{
    // The parameters 'body', 'location', 'otherBody' are all unused.
    // Return true if Moon is at greatest distance from Earth.
    
    var day1b = day1 + 1.0e-7;
    var day2a = day2 - 1.0e-7;

    var dist1  = Astronomy.Moon.DistanceFromEarth(day1);
    var dist1b = Astronomy.Moon.DistanceFromEarth(day1b);
    var dist2a = Astronomy.Moon.DistanceFromEarth(day2a);
    var dist2  = Astronomy.Moon.DistanceFromEarth(day2);

    return (dist1b <= dist1) && (dist2a <= dist2);
}

function Astronomy_MinDistance (body, day1, day2, location, otherBody)
{
    var day1b = day1 + 1.0e-7;
    var day2a = day2 - 1.0e-7;
}


function Astronomy_MinAngleWithOtherBodyCondition (body, day1, day2, location, otherBody)
{
    var day1b = day1 + 1.0e-7;
    var day2a = day2 - 1.0e-7;
    
    var sa1  = Astronomy.AngleBetweenBodiesInDegrees (body, otherBody, day1 );
    var sa1b = Astronomy.AngleBetweenBodiesInDegrees (body, otherBody, day1b);
    var sa2a = Astronomy.AngleBetweenBodiesInDegrees (body, otherBody, day2a);
    var sa2  = Astronomy.AngleBetweenBodiesInDegrees (body, otherBody, day2 );
    
    return (sa1 >= sa1b) && (sa2a <= sa2);
}

function Astronomy_PeakVisualMagnitudeCondition (body, day1, day2)
{
    var day1b = day1 + 1.0e-7;
    var day2a = day2 - 1.0e-7;
    
    var m1  = body.VisualMagnitude (day1);
    var m1b = body.VisualMagnitude (day1b);
    var m2a = body.VisualMagnitude (day2a);
    var m2  = body.VisualMagnitude (day2);
    
    return (m1 >= m1b) && (m2a <= m2);
}

function Astronomy_VernalEquinoxCondition (sun, day1, day2)
{
    var gc1 = sun.GeocentricCoordinates (day1);
    var gc2 = sun.GeocentricCoordinates (day2);
    return gc1.x > 0 && gc2.x > 0 && gc1.y <= 0 && gc2.y >= 0;
}

function Astronomy_AutumnalEquinoxCondition (sun, day1, day2)
{
    var gc1 = sun.GeocentricCoordinates (day1);
    var gc2 = sun.GeocentricCoordinates (day2);
    return gc1.x < 0 && gc2.x < 0 && gc1.y >= 0 && gc2.y <= 0;
}

function Astronomy_NorthernSolsticeCondition (sun, day1, day2)
{
    var gc1 = sun.GeocentricCoordinates (day1);
    var gc2 = sun.GeocentricCoordinates (day2);
    return gc1.y > 0 && gc2.y > 0 && gc1.x >= 0 && gc2.x <= 0;
}

function Astronomy_SouthernSolsticeCondition (sun, day1, day2)
{
    var gc1 = sun.GeocentricCoordinates (day1);
    var gc2 = sun.GeocentricCoordinates (day2);
    return gc1.y < 0 && gc2.y < 0 && gc1.x <= 0 && gc2.x >= 0;
}

function Astronomy_RelativeLongitudeCondition (body, day1, day2, location, longitude)
{
    var lon1 = Astronomy.RelativeLongitude (body, day1);
    var lon2 = Astronomy.RelativeLongitude (body, day2);
    return AnglesInOrder (lon1, longitude, lon2);
}


function AnglesInOrder (a, b, c)
{
    a = Angle.FixDegrees (a);
    b = Angle.FixDegrees (b);
    c = Angle.FixDegrees (c);

    var MARGIN = 45.0;
    if (c < a) {
        var swap = a;
        a = c;
        c = swap;
    }

    if (a <= MARGIN && c >= 360.0 - MARGIN) {
        // this is the wraparound case
        return ((b >= 0.0) && (b <= a)) || ((b >= c) && (b <= 360.0));
    } else if (c - a <= 2 * MARGIN) {
        return (a <= b) && (b <= c);
    } else {
        throw "AnglesInOrder continuity failure!";
    }
}



function Astronomy_FindNextTransition (body, startDay, endDay, location, condition, numIntervals, customData)
{
    var dt = (endDay - startDay) / numIntervals;
    var t1 = startDay;
    for (var i=0; i < numIntervals; ++i) {
        var t2 = t1 + dt;
        if (condition (body, t1, t2, location, customData)) {
            if (dt < 1.0e-6) {
                return t1 + (0.5 * dt);
            } else {
                var evt = Astronomy_FindNextTransition (body, t1, t2, location, condition, 3, customData);
                if (evt == null) {
                    return t1 + (0.5 * dt);    // somehow the answer slipped through the cracks: so do the best we can
                } else {
                    return evt;
                }
            }
        }
        
        t1 = t2;
    }
    return null;
}

//----------------------------------------------------------------------------------------------

function CartesianCoordinates (x, y, z)
{
    this.x = x;
    this.y = y;
    this.z = z;
}

CartesianCoordinates.prototype.toString = function()
{
    return "(" + this.x + ", " + this.y + ", " + this.z + ")";
}

CartesianCoordinates.prototype.Distance = function()
{
    return Math.sqrt (this.x*this.x + this.y*this.y + this.z*this.z);
}

CartesianCoordinates.prototype.Normalize = function()
{
    var r = this.Distance();
    if (r == 0.0) {
        throw "Cannot normalize zero vector.";
    }
    return new CartesianCoordinates (this.x/r, this.y/r, this.z/r);
}

CartesianCoordinates.prototype.Subtract = function (other)
{
    return new CartesianCoordinates (this.x - other.x, this.y - other.y, this.z - other.z);
}

CartesianCoordinates.prototype.Add = function (other)
{
    return new CartesianCoordinates (this.x + other.x, this.y + other.y, this.z + other.z);
}

function AngleBetweenVectorsInDegrees (va, vb)
{
    var a = va.Normalize ();
    var b = vb.Normalize ();

    // The normalized dot product of two vectors is equal to the cosine of the angle between them.
    var dotprod = a.x*b.x + a.y*b.y + a.z*b.z;
    var abs = Math.abs (dotprod);
    if (abs > 1.0) {
        if (abs > 1.00000001) {
            // Something is really messed up!  Avoid impossible Math.Acos() argument...
            throw "Internal error: dot product = " + dotprod;
        } else {
            // Assume that roundoff error has caused a valid dot product value to wander outside allowed values for Math.Acos()...
            return (dotprod < 0) ? 180.0 : 0.0;
        }
    } else {
        return Angle.DEG_FROM_RAD * Math.acos (dotprod);
    }
}

//----------------------------------------------------------------------------------------------

function SphericalCoordinates (longitude, latitude, radius)
{
    this.longitude = longitude;
    this.latitude  = latitude;
    this.radius    = radius;
}

//----------------------------------------------------------------------------------------------

function GeographicCoordinates (longitude, latitude, elevationInMetersAboveSeaLevel)
{
    if (elevationInMetersAboveSeaLevel == null) {
        elevationInMetersAboveSeaLevel = 0.0;
    }

    this.longitude = longitude;
    this.latitude  = latitude;
    this.radius    = (elevationInMetersAboveSeaLevel / METERS_PER_ASTRONOMICAL_UNIT) + EARTH_RADII_PER_ASTRONOMICAL_UNIT;
}

//----------------------------------------------------------------------------------------------
// Create a singleton instance "Astronomy", as a substitute for a namespace...

var Astronomy = new AstronomyClass();

module.exports = { Astronomy, Angle, GeographicCoordinates }

/*
    $Log: astronomy.js,v $
    Revision 1.33  2009/03/08 00:02:10  Don.Cross
    My C# astro.exe (compiled as sun.exe) now has a "javascript" option that generates constellation.js.
    This new constellation.js file contains the data needed for constellation calculation based on equatorial coordinates.
    I updated my astronomy.js code to allow use of this data to determine a constellation.
    The page solar_system.html uses this now to show the concise constellation symbol for each celestial body.

    Revision 1.32  2008/03/19 21:43:36  Don.Cross
    Added ability to calculate date/time of peak visual magnitude.
    Fixed bug where script crashed in December months due to borked up line of code.

    Revision 1.31  2008/03/19 20:43:55  Don.Cross
    Implemented visual magnitude calculation.
    Solar System page displays both magnitude and sun angle now.

    Revision 1.30  2008/03/17 18:27:40  Don.Cross
    Calculate conjunctions of planets.  For inferior planets, distinguish between inferior and superior conjunctions.

    Revision 1.29  2008/03/17 01:05:41  Don.Cross
    Implemented max elongation formulas for Mercury and Venus.
    I had to tweak the Astronomy_FindNextTransition to return a less precise answer rather than just give up
    when the answer slips through the cracks due to roundoff error.  This seems to work fine, but I will need to keep an eye on it.

    Revision 1.28  2008/03/16 22:54:47  Don.Cross
    Moved next/prev/current controls above the calendar so they don't jump around.
    Now calculate oppositions of superior planets.

    Revision 1.27  2008/03/16 21:59:37  Don.Cross
    Implemented calculation of equinoxes and solstices.
    Removed fixed height of extra-info div so that it expands when there are a lot of events on one day.
    Added links to other astronomy pages.

    Revision 1.26  2008/03/16 20:40:09  Don.Cross
    Now calculating moon phases consistently with SunCalc.exe (and Sky & Telescope magazine web site)
    using relative longitude with Sun instead of phase angles.
    Displaying full moon and new moon also.

    Revision 1.25  2008/03/16 20:07:07  Don.Cross
    Implemented calculation of moon's first and third quarters.
    Astronomy Calendar now displays graphics and adds details in detail box.

    Revision 1.24  2008/03/01 15:19:58  Don.Cross
    Added Astronomy.IsBodyName().
    Instead of ugly links from planet names to calculate extra info, made select box
    in the extra info box itself to choose the object.
    Made the choice persistent via a cookie.

    Revision 1.23  2008/02/28 23:28:11  Don.Cross
    Implemented calculation of rise, set, and culmination of a body.
    When user clicks on a body's name, display rise, set and culmination times for it.
    Default to showing sunrise, sun culmination, and sunset.

    Revision 1.22  2008/02/27 21:36:12  Don.Cross
    I have changed my mind about correcting altitude based on the apparent angular radius of the Sun and Moon.
    It makes more sense to calculate apparent altitude of the center of all objects, but then to know that the Sun and Moon
    have a significant angular radius for the purposes of rise and set times.
    This also allows the angular size of the Sun and Moon to be calculated exactly based on their distance from the observer.
    Also, because Civil, Nautical, and Astronomical Twilight are defined by how far below the horizon the center of the Sun is,
    I constricted the atmospheric correction to 6 degrees instead of 10 degrees.  Civil Twilight is when the Sun is at -6 degrees altitude,
    so I will be able to calculate that in agreement with other sources.

    Revision 1.21  2008/02/27 02:27:49  Don.Cross
    Oops!  Forgot to give a horizon correction for Pluto!
    This fixes a problem where altitude of Pluto shows up as "NaN" when it was within 10 degrees of the horizon.

    Revision 1.20  2008/02/26 23:17:20  Don.Cross
    Implemented refinements to topocentric parallax correction, factoring in the oblateness of the Earth.

    Revision 1.19  2008/02/22 21:52:18  Don.Cross
    Implemented Pluto.

    Revision 1.18  2008/02/21 12:40:35  Don.Cross
    No code change... just made indentation consistent on 1 line of code.

    Revision 1.17  2008/02/18 00:05:35  Don.Cross
    Fixed division by zero in Horizontal coordinate conversion when latitude = 0.

    Revision 1.16  2008/02/17 23:19:22  Don.Cross
    Got horizontal coordinates working!
    Fixed bugs in solar_system.html with not loading/applying polarity of lat/lon correctly.

    Revision 1.15  2008/02/17 20:07:16  Don.Cross
    Starting to add options to display different kinds of angular coordinates.
    Adding controls to enter geographic coordinates.

    Revision 1.14  2008/02/16 02:59:59  Don.Cross
    Display RA and DEC in nice HTML format with hours/degrees, minutes, seconds.

    Revision 1.13  2008/02/16 02:33:56  Don.Cross
    Implemented calculation of equatorial coordinates (RA, DEC), but not yet with topocentric parallax correction.

    Revision 1.12  2008/02/15 22:39:15  Don.Cross
    Renamed "Astronomy.Planet" array to "Astronomy.Body", and also includes in it the Sun and Moon.
    In the future I may include some famous comets and large asteroids also.

    Revision 1.11  2008/02/15 21:14:10  Don.Cross
    Implemented the Moon.
    Implemented geocentric coordinates for all bodies.

    Revision 1.10  2008/02/15 20:25:23  Don.Cross
    Added Sun as a celestial body.

    Revision 1.9  2008/02/15 16:56:13  Don.Cross
    Fixed quirk in Astronomy.DayValue function: passing 0 as the day value
    to Date.UTC is what was causing inconsistency in cscript and html version
    of Javascript.  When I pass 1 instead, they both return the same, correct value.

    Revision 1.8  2008/02/15 15:03:20  Don.Cross
    Implemented Neptune.
    Implemented perturbation algorithms, which allowed implementing Jupiter, Saturn, and Uranus.

    Revision 1.7  2008/02/15 03:01:25  Don.Cross
    Implemented Mercury and Venus.
    Display coordinates in a table instead of clunky text.
    Display coordinates updated every second in real time.

    Revision 1.6  2008/02/14 23:25:56  Don.Cross
    Fixed bug in AstronomyClass: was off by 1 day.
    Got Mars working, which is the framework for Mercury and Venus as well.
    To do the outer planets (except Pluto), I will need to implement perturbation.
    Pluto is a special case and will be coded separately.

    Revision 1.5  2008/02/14 22:36:05  Don.Cross
    Display human-friendly date and time along with numeric day value.

    Revision 1.4  2008/02/14 22:28:11  Don.Cross
    Got Earth cartesian coordinate formulas working.

    Revision 1.3  2008/02/14 22:04:43  Don.Cross
    Defined CartesianCoordinates class.
    Discovered it is easy to override toString(), which in turn overrides default behavior for converting to string.

    Revision 1.2  2008/02/14 21:43:40  Don.Cross
    Implemented Astronomy.DayValue() and Astronomy.CurrentDayValue().

    Revision 1.1  2008/02/14 21:20:32  Don.Cross
    I figured out how to do the equivalent of #include from cscript.

*/
