/*************************************************************************
 *     Peaks from array
 *     Copyright (C) agosto 2022 Esteban del Castillo Pérez
 *     esteban.delcastillo@urv.cat
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/

#include "peaks.h"

//Constructor
//arguments:
//magnitude_p[] -> array of input values to extract peaks (positive reals)
//magnitudeSize -> size of the array of input values
//peaksGap -> minimum increment between consecutive scans so that a peak can be considered as such
//noise -> minimum value so that the value of magnitude_p[] can be considered
Peaks::Peaks(float *magnitude_p, int magnitudeSize, float resolution, float maxSensitivity, float noise)
    {
    m_magnitude_p=magnitude_p;
    m_magnitudeSize=magnitudeSize;
    m_minGap=resolution;        //variable depending on magnitude.
    m_maxGap=maxSensitivity;
    m_noise=noise;
    m_mzIndex_p=0;
    m_mzIndex_p=new ION_INDEX[magnitudeSize];
    }

//destructor
//reserved memory is released.
Peaks::~Peaks()
    {
    if(m_mzIndex_p) delete [] m_mzIndex_p;
    }

//update the list of magnitudes
//arguments:
//  magnitude_p[] -> array of input values to extract peaks (positive reals).
//  magnitudeSize -> size of the array of input values.
void Peaks::setMagnitude(float *magnitude_p, int magnitudeSize)
{
    m_magnitude_p=magnitude_p;
    m_magnitudeSize=magnitudeSize;
}

//the sets of mz that make up an ion (mzLow, mzHigh and mzMax) are extracted in m_mzIndex_p.
//It does this considering two parameters: the noise at the input and a minimum increase with respect to previous values so that it can be considered a peak.
//After having the list, the limits of each ion are adjusted based on the noise.
//return the number of ions
//
//Mealy state machine (FSM)
//4 states: UP, DOWN, T_UP and T_DOWN
//UP status
//remains when positive increments greater than or equal to m_peaksGap occur.
//if positive increments less than m_peaksGap occur, it goes to the T_UP state.
//if negative increments greater than m_peaksGap occur, it is passed to T_DOWN (peak).
//DOWN status
//remains when negative increments greater than or equal to m_peaksGap occur.
//if negative increments less than m_peaksGap occur, the state goes to T_DOWN.
//if positive increments greater than m_peaksGap occur, it is passed to the UP (valley).
//T_UP status.
//remains when absolute increments occur with respect to the start less than or equal to m_peaksGap.
//If these increments are positive and exceed m_peaksGap, it goes to the UP state.
//If these increments are negative and exceed m_peaksGap, it goes to the DOWN state (peak at maximum value within the state).
//T_DOWN state.
//remains when absolute increments occur with respect to the start less than or equal to m_peaksGap.
//If these increments are positive and exceed m_peaksGap, it goes to the UP state (valley at minimum value within the state).
//If these increments are negative and exceed m_peaksGap, it goes to the DOWN state.

//Grades:
//when leaving the T-UP or TDOWN states, the last two elements must be reevaluated in the next state.
//The FSM is entered through the UP or DOWN states.
//When reaching the last element, special consideration is made about the state it is in to complete the last ion.
//Those elements with a value less than m_noise are canceled before analyzing them.
//Arguments
//mzIndexIni -> low index to the magnitude_p[] array. First element to evaluate.
//mzIndexEnd -> high index to the magnitude_p[] array. Last element to evaluate.
int Peaks::get(int mzIndexIni, int mzIndexEnd)
    {
    if(mzIndexEnd==0 || mzIndexEnd>=m_magnitudeSize) {return -1;}
    int mzSize=mzIndexEnd-mzIndexIni+1;

    bool hitCresta=false, hitValle=false;
    State state;
    int ionCresta=0, ionValle=0, ionIndex=0, TmaxIon=0, TminIon=0;
    float TiniVal=0, TmaxVal=0, TminVal=0, A, B, resolution;
;
    float delta;
    int indexIni=0;
//    Common tools;

    if(mzSize==1){  //case there is only one value.
            if(m_magnitude_p[0]<m_noise) return 0; //no peak detected.
            m_mzIndex_p[ionIndex].low= mzIndexIni;
            m_mzIndex_p[ionIndex].high=mzIndexIni;
            m_mzIndex_p[ionIndex].max=mzIndexIni;
            return 1; //un pico
            }
    else  //rest of cases (more than one value).
        {
        //the state of entry to the machine is determined.
        for(indexIni=0; indexIni<=mzIndexEnd; indexIni++)
            {
            A=m_magnitude_p[indexIni]; //noise is considered by canceling the lower values.
            if(A<m_noise) A=0.0;
            B=m_magnitude_p[indexIni+1];
            if(B<m_noise) B=0.0;

            //updates the resolution according to the level of information.
            resolution=rangeConversion(A>B?A:B, 0.0, 1.0, m_minGap, m_maxGap);
            //initial state for FSM.
            delta=B-A;
            if(delta>resolution)  {state=UP;   break;}
            else if(-delta>resolution) {state=DOWN; ionCresta=indexIni; hitCresta=true; break;}
            }
        if(indexIni>=mzIndexEnd) return 0; //there are no values to consider.

        //FSM operation.
        bool ionOK;
        for(int mzIndex=indexIni; mzIndex<=mzIndexEnd; mzIndex++)
            {
            ionOK=false;
            A=m_magnitude_p[mzIndex]; //noise is considered by canceling the lower values.
            if(A<m_noise) A=0.0;
            B=m_magnitude_p[mzIndex+1];
            if(B<m_noise) B=0.0;

            //updates the resolution according to the level of information.
            resolution=rangeConversion(A>B?A:B, 0.0, 1.0, m_minGap, m_maxGap);
            delta=B-A;
            if(state==UP)
                {
                if(delta>resolution) continue; //keep going up.
                else if(-delta>resolution) {state=DOWN; ionCresta=mzIndex; hitCresta=true;} //subía y ahora baja
                else  {state=T_UP; TiniVal=A; TmaxVal=TiniVal; TmaxIon=mzIndex;}
                }
            else if(state==DOWN)
                {
                if(-delta>resolution) continue; //it keeps going down.
                else if(delta>resolution) {state=UP; ionValle=mzIndex; hitValle=true;}
//                else if(delta<=m_noise)   {state=T_DOWN; ionValle=mzIndex; hitValle=true;} //finished peak to reach the noise level.
                else  {state=T_DOWN; TiniVal=A; TminVal=TiniVal; TminIon=mzIndex;}
                }
            else if(state==T_UP)
                {
                delta=A-TiniVal;
                if(A>=TmaxVal) {TmaxVal=A; TmaxIon=mzIndex;}

                if(delta>resolution) {state=UP; mzIndex--;} //comes out of the tube rising.
                else if(-delta>resolution) {state=DOWN; ionCresta=TmaxIon; hitCresta=true; mzIndex--;} //comes out of the tube going down.
                }
            else if(state==T_DOWN)
                {
                delta=A-TiniVal;
                if(A<=TminVal) {TminVal=A; TminIon=mzIndex;}

                if(delta>resolution) {state=UP; ionValle=TminIon; hitValle=true; mzIndex--;} //comes out of the tube rising.
                else if(-delta>resolution) {state=DOWN; mzIndex--;} //comes out of the tube going down.
                }
            //An ion is complete if there is a peak and valley on the right,
            //the valley on the left coincides with the valley on the right of the previous ion.
            if(hitValle && hitCresta)
                {
                if(ionIndex==0)
                    m_mzIndex_p[ionIndex].low=indexIni;
                else
                    m_mzIndex_p[ionIndex].low=m_mzIndex_p[ionIndex-1].high;
                m_mzIndex_p[ionIndex].high=ionValle;
                m_mzIndex_p[ionIndex].max=ionCresta;
                hitCresta=false; hitValle=false;
                state=UP;
                ionOK=true;
                ionIndex++;
                }
            }
        //end FSM
        //Right when it comes out the ion may not be well finished (without a valley). In this case, the spike may or may not be considered.
        //It is considered if it is down, there is a peak and its final magnitude is less than or equal to 1/4 of the maximum
        if(!ionOK)
            {
            if(state==DOWN && hitCresta && m_magnitude_p[mzIndexEnd]<=m_magnitude_p[ionCresta]/4)
                {
                if(ionIndex) //hay más ionIndexes antes
                    m_mzIndex_p[ionIndex].low=m_mzIndex_p[ionIndex-1].high;
                else
                    m_mzIndex_p[ionIndex].low=indexIni;
                m_mzIndex_p[ionIndex].high=mzIndexEnd;
                m_mzIndex_p[ionIndex].max=ionCresta;
                ionIndex++;
                }
            else if(state==T_DOWN && hitCresta && m_magnitude_p[mzIndexEnd]<=m_magnitude_p[ionCresta]/4)
                {
                if(ionIndex) //hay más ionIndexes antes
                    m_mzIndex_p[ionIndex].low=m_mzIndex_p[ionIndex-1].high;
                else
                    m_mzIndex_p[ionIndex].low=indexIni;
                m_mzIndex_p[ionIndex].high=mzIndexEnd;
                m_mzIndex_p[ionIndex].max=ionCresta;
                ionIndex++;
                }
            }
        }
     //up to here, an ion extends from the end of the previous one to the valley on the right of the next one.
     //Low and high limits are now adjusted when considering noise: an ion is delimited by two nearby valleys
     for(int i=0; i<ionIndex; i++) //for all ions.
        {
        int A=-1, B=-1;
        int max=m_mzIndex_p[i].max;
        for(int j=m_mzIndex_p[i].low; j<m_mzIndex_p[i].high; j++) //for the values within the ion.
            {
            if(j<max && m_magnitude_p[j]<=m_noise) A=j; //before the maximum.
            else if(j>max && m_magnitude_p[j]<=m_noise) {B=j; break;}  //after the maximum.
            }
        if(A!=-1) m_mzIndex_p[i].low=A; //adjustment.
        if(B!=-1) m_mzIndex_p[i].high=B;
        }
    //We look to see if the magnitude differences between valleys exceed resolution.
    for(int i=0; i<ionIndex; i++) //for all ions.
        {
        if(m_mzIndex_p[i].high - m_mzIndex_p[i].low <= 2) //one peak and two valleys or less.
            {m_mzIndex_p[i].confidence=true; continue;}

        m_mzIndex_p[i].confidence=false;
        float minValue=1e32, maxValue=-1;
        for(int j=m_mzIndex_p[i].low+1; j<m_mzIndex_p[i].high; j++) //for the values within the ion.
            {
            if(m_magnitude_p[j]>maxValue) maxValue=m_magnitude_p[j];
            if(m_magnitude_p[j]<minValue) minValue=m_magnitude_p[j];
            }
        //updates the resolution according to the level of information.
        resolution=rangeConversion(maxValue, 0.0, 1.0, m_minGap, m_maxGap);
        if(maxValue-minValue>resolution)
                m_mzIndex_p[i].confidence=true;
        }
    return ionIndex;
}

//Applies a linear conversion between ranges
//given the value 'value' located within the range [x1, x2] returns its equivalent within
//of the range [y1, y2]. The points (x1, y1), (x2, y2) are known
//Saturates if 'value' is out of range or if the target range is null
double Peaks::rangeConversion(double value, double x1, double x2, double y1, double y2)
{
  if(value<x1)
      return y1;
  else if(value>x2)
      return y2;
  else if(y2-y1==0)
      return y1;

  return y1+value*(y2-y1)/(x2-x1);
}
