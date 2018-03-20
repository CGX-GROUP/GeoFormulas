/** \file CrsIntersect.cpp
*   \brief
*/

/****************************************************************************/
/*  CrsIntersect.cpp                                                        */
/****************************************************************************/
/*                                                                          */
/*  Copyright 2008 - 2016 Paul Kohut                                        */
/*  Licensed under the Apache License, Version 2.0 (the "License"); you may */
/*  not use this file except in compliance with the License. You may obtain */
/*  a copy of the License at                                                */
/*                                                                          */
/*  http://www.apache.org/licenses/LICENSE-2.0                              */
/*                                                                          */
/*  Unless required by applicable law or agreed to in writing, software     */
/*  distributed under the License is distributed on an "AS IS" BASIS,       */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         */
/*  implied. See the License for the specific language governing            */
/*  permissions and limitations under the License.                          */
/*                                                                          */
/****************************************************************************/

#include "Conversions.h"
#include "GeoFormulas.h"


namespace GeoCalcs {
    /**
    *
    */
    bool CrsIntersect(const LLPoint &llPt1, double az13,
                      const LLPoint &llPt2, double az23, double dTol, LLPoint &llIntersect)
    {
        double az31, dist13, az32, dist23;
        az31 = az32 = 0.0;
        dist13 = dist23 = 0.0;
        return CrsIntersect(llPt1, az13, az31, dist13,
                            llPt2, az23, az32, dist23,
                            dTol, llIntersect);

    }


    /**
    *
    */
    bool CrsIntersect(const LLPoint &llPt1, double az13,
                      double &az31, double &dist13, const LLPoint &llPt2, double az23,
                      double &az32, double &dist23, double dTol, LLPoint &llIntersect)
    {
        LLPoint pt1 = llPt1;
        LLPoint pt2 = llPt2;
        double dAz13 = az13;
		double dAz23 = az23;
        InverseResult result;
        DistVincenty(pt1, pt2, result);
        double dist12 = result.distance;
		double crs12 = result.azimuth;
		double crs21 = result.reverseAzimuth;

        double angle1 = SignAzimuthDifference(crs12, dAz13);
		double angle2 = SignAzimuthDifference(dAz23, crs21);

        if (sin(angle1) == 0.0 && sin(angle2) == 0.0)
            return false;

        if(sin(angle1)*sin(angle2)<0)
		{
			//The courses lay on opposite sides of the line between p1 and p2. One of the course must be reversed.
			if(fabs(angle1)>fabs(angle2))
			{
				dAz13 = Mod(dAz13 + M_PI, M_2PI);
				angle1 = SignAzimuthDifference(crs12, dAz13);
			}
			else
			{
				dAz23 = Mod(dAz23 + M_PI, M_2PI);
				angle2 = SignAzimuthDifference(dAz23, crs21);
			}
		}
		
		angle1 = fabs(angle1);
		angle2 = fabs(angle2);

        double cosA = cos(angle1);
        double sinA = sin(angle1);
        double cosB = cos(angle2);
        double sinB = sin(angle2);
        
        //Test if the intersection is backcourse. In this case, reverse both courses.
		//We check the sign of the Napier's analogies result (cf. https://en.wikipedia.org/wiki/Solution_of_triangles#A_side_and_two_adjacent_angles_given_.28spherical_ASA.29 )
		double tanc2 = tan(dist12 / (2*kSphereRadius));
		double cotc2 = 1/tanc2;
		double sinAplusB = sin(angle1+angle2);
		double sinAminusB = sin(angle1-angle2);
        
		double b_ = kSphereRadius * atan( 2*sinB/(cotc2*sinAplusB+tanc2*sinAminusB) );
		if(sinB/(cotc2*sinAplusB+tanc2*sinAminusB) <= 0) //ie. b <= 0
		{
			dAz13 = Mod(dAz13 + M_PI, M_2PI);
			angle1 = fabs(SignAzimuthDifference(crs12, dAz13));
			dAz23 = Mod(dAz23 + M_PI, M_2PI);
			angle2 = fabs(SignAzimuthDifference(dAz23, crs21));
		
			cosA = cos(angle1);
			sinA = sin(angle1);
			cosB = cos(angle2);
			sinB = sin(angle2);
		}

		// step 7
		// locate approx intersection of point3 using spherical earth model. Section 3.2

        const double C = acos(-cosA * cosB + sinA * sinB * cos(dist12 / kSphereRadius));
        const double cosC = cos(C);
        const double sinC = sin(C);
        const double a = kSphereRadius * acos((cosA + cosB * cosC) / (sinB * sinC));
        const double b = kSphereRadius * acos((cosB + cosA * cosC) / (sinA * sinC));

        if (std::isnan(a) || std::isnan(b))
            return false;

        llIntersect = DestVincenty(pt1, dAz13, b);
        DistVincenty(pt1, llIntersect, result);
        dist13 = result.distance;

        LLPoint llInv = llIntersect;
        llInv.latitude = -llInv.latitude;
        llInv.longitude = llInv.longitude + M_PI - M_2PI;
        DistVincenty(pt1, llInv, result);

        if (dist13 > result.distance)
        {
            llIntersect = llInv;
            dist13 = result.distance;
            az31 = result.reverseAzimuth;
            dAz13 = dAz13 + M_PI;
            dAz23 = dAz23 + M_PI;
        }

        DistVincenty(pt2, llIntersect, result);
        dist23 = result.distance;

        if (dist13 < NmToMeters(1))
        {
            pt1 = DestVincenty(pt1, dAz13 + M_PI, NmToMeters(1.0));
            DistVincenty(pt1, llIntersect, result);
            dAz13 = result.azimuth;
        }
        if (dist23 < NmToMeters(1))
        {
            pt2 = DestVincenty(pt2, dAz23 + M_PI, NmToMeters(1.0));
            DistVincenty(pt2, llIntersect, result);
            dAz23 = result.azimuth;
        }

        bool bSwapped = false;
        if (dist23 < dist13)
        {
            std::swap(pt1, pt2);
            std::swap(dAz13, dAz23);
            dist13 = dist23;
            bSwapped = true;
        }

        double distarray[2], errarray[2];
        distarray[0] = dist13;
        distarray[1] = 1.01 * dist13;

        llIntersect = DestVincenty(pt1, dAz13, dist13);
        DistVincenty(pt2, llIntersect, result);
        errarray[0] = SignAzimuthDifference(result.azimuth, dAz23);

        llIntersect = DestVincenty(pt1, dAz13, distarray[1]);
        DistVincenty(pt2, llIntersect, result);
        errarray[1] = SignAzimuthDifference(result.azimuth, dAz23);

        const int nMaxCount = 15;
        double dErr = 0.0;
        int k = 0;
        while (k == 0 || ((dErr > dTol) && (k <= nMaxCount)))
        {
            FindLinearRoot(distarray, errarray, dist13);
            if (std::isnan(dist13))
                break;

            LLPoint newPt = DestVincenty(pt1, dAz13, dist13);
            DistVincenty(pt2, newPt, result);
            errarray[0] = errarray[1];
            errarray[1] = SignAzimuthDifference(result.azimuth, dAz23);
            distarray[0] = distarray[1];
            distarray[1] = dist13;

            DistVincenty(newPt, llIntersect, result);
            dErr = result.distance;
            llIntersect = newPt;
            k++;
        }

        // display if k == maxinteratorcount (10) and show error message
        // because results might not have converged.
        if (k > nMaxCount && dErr > 1e-8)
            return false;

        if (bSwapped)
        {
            std::swap(pt1, pt2);
            std::swap(dAz13, dAz23);
            dist13 = dist23;
        }
        DistVincenty(llIntersect, pt1, result);
        az31 = result.azimuth;
        dist13 = result.distance;

        DistVincenty(llIntersect, pt2, result);
        az32 = result.azimuth;
        dist23 = result.distance;

        return true;
    }
}
