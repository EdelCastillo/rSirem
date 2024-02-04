/*************************************************************************
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

#ifndef GET_PEAKS_H
#define GET_PEAKS_H


//Determines the existing peaks in an array of magnitudes whose values must be in the range [0:1]. 
//A peak is made up of three elements: valley on the left, crest and valley on the right. The sensitivity
//is variable within a given range and is directly associated with the value of the magnitude.
class Peaks
{
public:
    enum State
        {
        UP,     //going up.
        DOWN,   //going down.
        T_UP,   //tube that you enter uphill.
        T_DOWN, //tube into which you enter downhill.
        };
        
    typedef struct //indices to the magnitude_p[] array that delimit an ion.
        {
        int low, high, max; //indexes to elements left valley, right valley and maximum.
        bool confidence;    //true if the magnitude variation between valleys exceeds m_peaksGap.
        }ION_INDEX;

    
    //Constructor
    //arguments:
    //magnitude_p[] -> array of input values to extract peaks (positive reals)
    //magnitudeSize -> size of the array of input values
    //peaksGap -> minimum increment between consecutive scans so that a peak can be considered as such
    //noise -> minimum value so that the value of magnitude_p[] can be considered
    Peaks(float *magnitude_p, int magnitudeSize, float minSensitivity, float maxSensitivity, float noise);
    
    //destructor
    //reserved memory is released.
    ~Peaks();
    
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
    int get(int mzIndexIni, int mzIndexEnd);

    //update the list of magnitudes
    //arguments:
    //  magnitude_p[] -> array of input values to extract peaks (positive reals).
    //  magnitudeSize -> size of the array of input values.
    void setMagnitude(float *magnitude_p, int magnitudeSize);
    
    //Applies a linear conversion between ranges
    //given the value 'value' located within the range [x1, x2] returns its equivalent within
    //of the range [y1, y2]. The points (x1, y1), (x2, y2) are known
    //Saturates if 'value' is out of range or if the target range is null
    double rangeConversion(double value, double x1, double x2, double y1, double y2);

    float   *m_magnitude_p, //array of input values to extract peaks.
            m_noise;        //minimum value so that the value of magnitude_p[] can be considered.
    int     m_magnitudeSize;//size of the array of input values.
    ION_INDEX *m_mzIndex_p; //ION_INDEX structure array where the peaks remain.
    float   m_minGap,       //minimum increment of magnitude between scans. Lower values are rejected (min sensitivity)
            m_maxGap;       //maximum increment of magnitude between scans. Higher values do not alter the result (max sensitivity)
};

#endif


