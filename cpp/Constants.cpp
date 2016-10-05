/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "Constants.h"
#include <math.h>

const double Constants::PI = 3.1415926535897932384626433832795;
const double Constants::Degree2Radian = PI / 180.0;
const double Constants::Arcmin2Radian = Degree2Radian / 60.0;
const double Constants::Arcsec2Radian = Arcmin2Radian / 60.0;
const double Constants::Radian2Degree = 1.0 / Degree2Radian;
const double Constants::Radian2Arcmin = 1.0 / Arcmin2Radian;
const double Constants::Radian2ArcSecond = 1.0 / Arcsec2Radian;
const double Constants::SquareRadian2SquareDegree = Radian2Degree * Radian2Degree;
const double Constants::WholeSphereInSquareDegree = 4.0 * PI * SquareRadian2SquareDegree;	
const double Constants::WholeSphereInSquareArcsec = 4.0 * PI * pow(180.0 / PI * 3600.0, 2);
