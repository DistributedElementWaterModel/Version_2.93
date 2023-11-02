/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments,          *
 *  river networks, lakes and landscape elements.                                       *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
enum MODEL_STRUCTURE { HBV_MODEL , KWA_MODEL };
enum LANDSURFACE { OPEN, BOG, FOREST, ALPINE, HEATHER, ROCK, GLACIER };
enum SOIL { OPEN_SOIL, PEAT, FOREST_SOIL, ALPINE_SOIL, HEATHER_SOIL, BEDROCK, GLACIER_BED };
//enum LANDSURFACE { OPEN, FOREST, ALPINE, HEATHER, ROCK, GLACIER };
//enum SOIL { OPEN_SOIL, FOREST_SOIL, ALPINE_SOIL, HEATHER_SOIL, BEDROCK, GLACIER_BED };
#define ELEMENT(a,b) (((a)*nCols)+(b))
const int numberModelStructures=2;           // All possible model structures
const int maximumNumberLandClasses=2;        // Maximum number of land/soil classes in use for each computational element
const int numberLandSurfaceClasses=7;        // All possible land surface types, excluding lakes
const int numberSoilClasses=7;               // All possible soil/subsurface types, excluding lakes
const int maximumCorrectionCatchments=1000;
const int numberInputSeries=2;
const int numberSnowClasses=9;
const int numberCharacteristic=100;          // Number of characteristic curves in kinematic wave model of saturated subsurface flow
const int minimumTimeStep=3600;
const int finalYear=2100;
const int finalMonth=12;
const int finalDay=31;
const int finalHour=23;
const int finalMinute=59;
const double missingData=-9999.0;
const double probNorm[9] = {0.01,0.04,0.1,0.2,0.3,0.2,0.1,0.04,0.01};
const double epsilon=1.0e-4;

class DistributedElement;
class WaterCourse;


//class ParametersCommon
class ParametersCommon
{
 public:
  ParametersCommon();
  ~ParametersCommon();
  void SetSECONDS_TIMESTEP(int value) { SECONDS_TIMESTEP = value; }
  int GetSECONDS_TIMESTEP() const { return SECONDS_TIMESTEP; }
  void SetNUM_PREC_SERIES(int value) { NUM_PREC_SERIES = value; }
  int GetNUM_PREC_SERIES() const { return NUM_PREC_SERIES; }
  void SetNUM_TEMP_SERIES(int value) { NUM_TEMP_SERIES = value; }
  int GetNUM_TEMP_SERIES() const { return NUM_TEMP_SERIES; }
  void SetPREC_GRAD_LOW(double value) { PREC_GRAD_LOW = value; }
  double GetPREC_GRAD_LOW() const { return PREC_GRAD_LOW; }
  void SetPREC_GRAD_HIGH(double value) { PREC_GRAD_HIGH = value; }
  double GetPREC_GRAD_HIGH() const { return PREC_GRAD_HIGH; }
  void SetGRAD_CHANGE_ALT(double value) { GRAD_CHANGE_ALT = value; }
  double GetGRAD_CHANGE_ALT() const { return GRAD_CHANGE_ALT; }
  void SetPREC_CORR_RAIN(double value) { PREC_CORR_RAIN = value; }
  double GetPREC_CORR_RAIN() const { return PREC_CORR_RAIN; }
  void SetPREC_CORR_SNOW(double value) { PREC_CORR_SNOW = value; }
  double GetPREC_CORR_SNOW() const { return PREC_CORR_SNOW; }
  void SetLAPSE_DRY(double value) { LAPSE_DRY = value; }
  double GetLAPSE_DRY() const { return LAPSE_DRY; }
  void SetLAPSE_WET(double value) { LAPSE_WET = value; }
  double GetLAPSE_WET() const { return LAPSE_WET; }
  void SetDAY_TEMP_MEMORY(double value) { DAY_TEMP_MEMORY = value; }
  double GetDAY_TEMP_MEMORY() const { return DAY_TEMP_MEMORY; }
  void SetLAKE_EPOT_PAR(double value) { LAKE_EPOT_PAR = value; }
  double GetLAKE_EPOT_PAR() const { return LAKE_EPOT_PAR; }
  void SetKLAKE(double value) { KLAKE = value; }
  double GetKLAKE() const { return KLAKE; }
  void SetDELTA_LEVEL(double value) { DELTA_LEVEL = value; }
  double GetDELTA_LEVEL() const { return DELTA_LEVEL; }
  void SetNLAKE(double value) { NLAKE = value; }
  double GetNLAKE() const { return NLAKE; }
  void SetMAXIMUM_LEVEL(double value) { MAXIMUM_LEVEL = value; }
  double GetMAXIMUM_LEVEL() const { return MAXIMUM_LEVEL; }
  void SetINITIAL_SOIL_MOISTURE(double value) { INITIAL_SOIL_MOISTURE = value; }
  double GetINITIAL_SOIL_MOISTURE() const { return INITIAL_SOIL_MOISTURE; }
  void SetINITIAL_UPPER_ZONE(double value) { INITIAL_UPPER_ZONE = value; }
  double GetINITIAL_UPPER_ZONE() const { return INITIAL_UPPER_ZONE; }
  void SetINITIAL_LOWER_ZONE(double value) { INITIAL_LOWER_ZONE = value; }
  double GetINITIAL_LOWER_ZONE() const { return INITIAL_LOWER_ZONE; }
  void SetINITIAL_SATURATED_ONE(double value) { INITIAL_SATURATED_ONE = value; }
  double GetINITIAL_SATURATED_ONE() const { return INITIAL_SATURATED_ONE; }
  void SetINITIAL_SATURATED_TWO(double value) { INITIAL_SATURATED_TWO = value; }
  double GetINITIAL_SATURATED_TWO() const { return INITIAL_SATURATED_TWO; }
  void SetINITIAL_LAKE_TEMP(double value) { INITIAL_LAKE_TEMP = value; }
  double GetINITIAL_LAKE_TEMP() const { return INITIAL_LAKE_TEMP; }
  void SetINITIAL_LAKE_LEVEL(double value) { INITIAL_LAKE_LEVEL = value; }
  double GetINITIAL_LAKE_LEVEL() const { return INITIAL_LAKE_LEVEL; }
  void SetINITIAL_SNOW(double value) { INITIAL_SNOW = value; }
  double GetINITIAL_SNOW() const { return INITIAL_SNOW; }
  void SetINITIAL_TOTAL_RESERVOIR(double value) { INITIAL_TOTAL_RESERVOIR = value; }
  double GetINITIAL_TOTAL_RESERVOIR() { return INITIAL_TOTAL_RESERVOIR; }
  void SetDAY_SNOW_ZERO(int value) { DAY_SNOW_ZERO = value; }
  int GetDAY_SNOW_ZERO() const { return DAY_SNOW_ZERO; }
  //  void SetNumStations(ifstream &fin, int numberPrecSeries, int numberTempSeries);
  //  double GetSTATION_ALTITUDE(int i) { return STATION_ALTITUDE[i]; }
  //  double GetSTATION_WEIGHT(int i) { return STATION_WEIGHT[i]; }
  
 private:
  int SECONDS_TIMESTEP;
  int NUM_PREC_SERIES;
  int NUM_TEMP_SERIES;
  int DAY_SNOW_ZERO;
  double PREC_GRAD_LOW;
  double PREC_GRAD_HIGH;
  double GRAD_CHANGE_ALT;
  double PREC_CORR_RAIN;
  double PREC_CORR_SNOW;
  double LAPSE_DRY;
  double LAPSE_WET;
  double DAY_TEMP_MEMORY;
  double LAKE_EPOT_PAR;
  double KLAKE;
  double DELTA_LEVEL;
  double NLAKE;
  double MAXIMUM_LEVEL;
  double INITIAL_SOIL_MOISTURE;
  double INITIAL_UPPER_ZONE;
  double INITIAL_LOWER_ZONE;
  double INITIAL_SATURATED_ONE;
  double INITIAL_SATURATED_TWO;
  double INITIAL_LAKE_TEMP;
  double INITIAL_LAKE_LEVEL;
  double INITIAL_SNOW;
  double INITIAL_TOTAL_RESERVOIR;
  //  double * STATION_WEIGHT;
  //  double * STATION_ALTITUDE;
};

ParametersCommon::ParametersCommon():
  SECONDS_TIMESTEP(0),
  NUM_PREC_SERIES(0),
  NUM_TEMP_SERIES(0),
  DAY_SNOW_ZERO(0),
  PREC_GRAD_LOW(0.0),
  PREC_GRAD_HIGH(0.0),
  GRAD_CHANGE_ALT(0.0),          
  PREC_CORR_RAIN(0.0),
  PREC_CORR_SNOW(0.0),
  LAPSE_DRY(0.0),
  LAPSE_WET(0.0),
  DAY_TEMP_MEMORY(0.0),
  LAKE_EPOT_PAR(0.0),
  KLAKE(0.0),
  DELTA_LEVEL(0.0),
  NLAKE(0.0),
  MAXIMUM_LEVEL(0.0),
  INITIAL_SOIL_MOISTURE(0.0),
  INITIAL_UPPER_ZONE(0.0),
  INITIAL_LOWER_ZONE(0.0),
  INITIAL_SATURATED_ONE(0.0),
  INITIAL_SATURATED_TWO(0.0),
  INITIAL_LAKE_TEMP(0.0),
  INITIAL_LAKE_LEVEL(0.0),
  INITIAL_SNOW(0.0),
  INITIAL_TOTAL_RESERVOIR(0.0)
{
}
     
ParametersCommon::~ParametersCommon()
{     
}

/*void ParametersCommon::SetNumStations(ifstream &fin, int numberPrecSeries, int numberTempSeries) 
  {
  int i;
  double sumWeight;
  STATION_WEIGHT = new double [numberPrecSeries+numberTempSeries];
  STATION_ALTITUDE = new double [numberPrecSeries+numberTempSeries];
  sumWeight=0.0;
  for (i=0; i<numberPrecSeries; i++) {
  fin.ignore(100,':'); 
  fin >> STATION_ALTITUDE[i];
  fin.ignore(100,':'); 
  fin >> STATION_WEIGHT[i];
  sumWeight = sumWeight+STATION_WEIGHT[i];
  }
  if (sumWeight != 1.0) {
  for (i=0; i<numberPrecSeries; i++) {
  STATION_WEIGHT[i] = STATION_WEIGHT[i]/sumWeight;
  }
  }
  sumWeight=0.0;
  for (i=0; i<numberTempSeries; i++) {
  fin.ignore(100,':'); 
  fin >> STATION_ALTITUDE[numberPrecSeries+i];
  fin.ignore(100,':'); 
  fin >> STATION_WEIGHT[numberPrecSeries+i];
  sumWeight = sumWeight+STATION_WEIGHT[numberPrecSeries+i];
  }
  if (sumWeight != 1.0) {
  for (i=0; i<numberTempSeries; i++) {
  STATION_WEIGHT[numberPrecSeries+i] = STATION_WEIGHT[numberPrecSeries+i]/sumWeight;
  }
  }
  }*/


//class ParametersLandSurface
class ParametersLandSurface
{
 public:
  ParametersLandSurface();
  ~ParametersLandSurface();
  void SetINTER_MAX(double value) { INTER_MAX = value; }
  double GetINTER_MAX() const { return INTER_MAX; }
  void SetEPOT_PAR(double value) { EPOT_PAR = value; }
  double GetEPOT_PAR() const { return EPOT_PAR; }
  void SetWET_PER_CORR(double value) { WET_PER_CORR = value; }
  double GetWET_PER_CORR() const { return WET_PER_CORR; }
  void SetACC_TEMP(double value) { ACC_TEMP = value; }
  double GetACC_TEMP() const { return ACC_TEMP; }
  void SetMELT_TEMP(double value) { MELT_TEMP = value; }
  double GetMELT_TEMP() const { return MELT_TEMP; }
  void SetSNOW_MELT_RATE(double value) { SNOW_MELT_RATE = value; }
  double GetSNOW_MELT_RATE() const { return SNOW_MELT_RATE; }
  void SetICE_MELT_RATE(double value) { ICE_MELT_RATE = value; }
  double GetICE_MELT_RATE() const { return ICE_MELT_RATE; }
  void SetFREEZE_EFF(double value) { FREEZE_EFF = value; }
  double GetFREEZE_EFF() const { return FREEZE_EFF; }
  void SetMAX_REL(double value) { MAX_REL = value; }
  double GetMAX_REL() const { return MAX_REL; }
  void SetALBEDO(double value) { ALBEDO = value; }
  double GetALBEDO() const { return ALBEDO; }
  void SetCV_SNOW(double value) { CV_SNOW = value; }
  double GetCV_SNOW() const { return CV_SNOW; }
  void SetSNOW_WEIGHT(int k, double value) { SNOW_WEIGHT[k]= value; }
  double GetSNOW_WEIGHT(int k) const { return SNOW_WEIGHT[k]; }
  
 private:
  double INTER_MAX;
  double EPOT_PAR;
  double WET_PER_CORR;
  double ACC_TEMP;
  double MELT_TEMP;  
  double SNOW_MELT_RATE;
  double ICE_MELT_RATE;
  double FREEZE_EFF;
  double MAX_REL;
  double ALBEDO;    
  double CV_SNOW;
  double SNOW_WEIGHT[numberSnowClasses];
};

ParametersLandSurface::ParametersLandSurface():
  INTER_MAX(0.0),
  EPOT_PAR(0.0),
  WET_PER_CORR(0.0),
  ACC_TEMP(0.0),        
  MELT_TEMP(0.0),       
  SNOW_MELT_RATE(0.0),  
  ICE_MELT_RATE(0.0),  
  FREEZE_EFF(0.0),      
  MAX_REL(0.0),         
  ALBEDO(0.0),
  CV_SNOW(0.0)
{
  int i;
  for (i=0; i<numberSnowClasses; i++) SNOW_WEIGHT[i]=0.0;
}
     
ParametersLandSurface::~ParametersLandSurface()
{     
}


//class ParametersSubSurfaceHbv
class ParametersSubSurfaceHbv
{
 public:
  ParametersSubSurfaceHbv();
  ~ParametersSubSurfaceHbv();
  void SetFC(double value) { FC = value; }
  double GetFC() const { return FC; }
  void SetFCDEL(double value) { FCDEL = value; }
  double GetFCDEL() const { return FCDEL; }
  void SetBETA(double value) { BETA = value; }
  double GetBETA() const { return BETA; }
  void SetINFMAX(double value) { INFMAX = value; }
  double GetINFMAX() const { return INFMAX; }
  void SetKUZ(double value) { KUZ = value; }
  double GetKUZ() const { return KUZ; }
  void SetALFA(double value) { ALFA = value; }
  double GetALFA() const { return ALFA; }
  void SetPERC(double value) { PERC = value; }
  double GetPERC() const { return PERC; }
  void SetKLZ(double value) { KLZ = value; }
  double GetKLZ() const { return KLZ; }
  void SetDRAW(double value) { DRAW = value; }
  double GetDRAW() const { return DRAW; }
  
 private:
  double FC;
  double FCDEL;
  double BETA;  
  double INFMAX;         
  double KUZ;         
  double ALFA;  
  double PERC;
  double KLZ;               
  double DRAW;              
};

ParametersSubSurfaceHbv::ParametersSubSurfaceHbv():
  FC(0.0),
  FCDEL(0.0),
  BETA(0.0),  
  INFMAX(0.0),         
  KUZ(0.0),         
  ALFA(0.0),  
  PERC(0.0),
  KLZ(0.0),               
  DRAW(0.0)
{
}
     
ParametersSubSurfaceHbv::~ParametersSubSurfaceHbv()
{     
}


//class ParametersKiWa
class ParametersKiWa
{
 public:
  ParametersKiWa();
  ~ParametersKiWa();
  //  void SetSLOPE_LENGTH(double value) { SLOPE_LENGTH = value; }
  //  double GetSLOPE_LENGTH() const { return SLOPE_LENGTH; }
  void SetSOIL_DEPTH(double value) { SOIL_DEPTH = value; }
  double GetSOIL_DEPTH() const { return SOIL_DEPTH; }
  void SetOV_PAR_1(double value) { OV_PAR_1 = value; }
  double GetOV_PAR_1() const { return OV_PAR_1; }
  void SetOV_PAR_2(double value) { OV_PAR_2 = value; }
  double GetOV_PAR_2() const { return OV_PAR_2; }
  void SetTSAT_0(double value) { TSAT_0 = value; }
  double GetTSAT_0() const { return TSAT_0; }
  void SetEFF_POR(double value) { EFF_POR = value; }
  double GetEFF_POR() const { return EFF_POR; }
  void SetKSAT_0(double value) { KSAT_0 = value; }
  double GetKSAT_0() const { return KSAT_0; }
  void SetA(double value) { A = value; }
  double GetA() const { return A; }
  void SetDELTA(double value) { DELTA = value; }
  double GetDELTA() const { return DELTA; }
  void SetLAMBDA_KW(double value) { LAMBDA_KW = value; }
  double GetLAMBDA_KW() const { return LAMBDA_KW; }
  void SetROOT_DEPTH(double value) { ROOT_DEPTH = value; }
  double GetROOT_DEPTH() const { return ROOT_DEPTH; }
  void SetWILT_POINT(double value) { WILT_POINT = value; }
  double GetWILT_POINT() const { return WILT_POINT; }
  void SetEACT_PAR(double value) { EACT_PAR = value; }
  double GetEACT_PAR() const { return EACT_PAR; }
  
 private:
  //  double SLOPE_LENGTH;
  double SOIL_DEPTH;  
  double OV_PAR_1;         
  double OV_PAR_2;         
  double TSAT_0;  
  double EFF_POR;         
  double KSAT_0;
  double A;               
  double DELTA;              
  double LAMBDA_KW;       
  double ROOT_DEPTH;      
  double WILT_POINT;           
  double EACT_PAR;
};

ParametersKiWa::ParametersKiWa():
  //  SLOPE_LENGTH(0),
  SOIL_DEPTH(0),  
  OV_PAR_1(0),         
  OV_PAR_2(0),         
  TSAT_0(0),  
  EFF_POR(0),         
  KSAT_0(0),
  A(0),               
  DELTA(0),              
  LAMBDA_KW(0),       
  ROOT_DEPTH(0),      
  WILT_POINT(0),           
  EACT_PAR(0.0)
{
}
     
ParametersKiWa::~ParametersKiWa()
{     
}


//class SelectedKiWaHillslopeElements
class SelectedKiWaHillslopeElements
{
 public:
  SelectedKiWaHillslopeElements();
  ~SelectedKiWaHillslopeElements();
  void SetNumberElements(int value) { numberElements = value; }
  int GetNumberElements() const { return numberElements; }
  void SetKiWaElement(int index, int value) { KiWaElements[index] = value; }
  int GetKiWaElement(int index) const { return KiWaElements[index]; }
  void SelectedKiWaHillslopeElementsInput(ifstream &fileControl, ofstream &fout);

 private:
  int numberElements;
  int * KiWaElements;
};

SelectedKiWaHillslopeElements::SelectedKiWaHillslopeElements():
  numberElements(0)
{
}
     
SelectedKiWaHillslopeElements::~SelectedKiWaHillslopeElements()
{     
}

void SelectedKiWaHillslopeElements::SelectedKiWaHillslopeElementsInput(ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  int i, numElements, elementId;

  /*  cout << " File with KiWa elements selected for hillslope output: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream fileKiWa(fileName);  // Open for reading
  if (fileKiWa == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  fileKiWa.ignore(100,':');
  fileKiWa >> numElements;
  SetNumberElements(numElements);
  KiWaElements = new int [numElements];
  for (i=0; i<numElements; i++) {
    fileKiWa >> elementId;
    SetKiWaElement(i, elementId);
    fileKiWa.ignore(256,'\n');
  }
  fout << endl << "Number of KiWa elements selected for hillslope output: " << GetNumberElements() << endl;
  for (i=0; i<GetNumberElements(); i++) {
    fout << GetKiWaElement(i) << endl;
  }
  fout << endl;
}


//class SelectedKiWaTimeSeriesElements
class SelectedKiWaTimeSeriesElements
{
 public:
  SelectedKiWaTimeSeriesElements();
  ~SelectedKiWaTimeSeriesElements();
  void SetNumberElements(int value) { numberElements = value; }
  int GetNumberElements() const { return numberElements; }
  void SetKiWaTimeSeriesElement(int index, int value) { timeSeriesElements[index] = value; }
  int GetKiWaTimeSeriesElement(int index) const { return timeSeriesElements[index]; }
  void SetLengthFractionOne(int index, double value) { lengthFractionOne[index] = value; }
  double GetLengthFractionOne(int index) const { return lengthFractionOne[index]; }
  void SetLengthFractionTwo(int index, double value) { lengthFractionTwo[index] = value; }
  double GetLengthFractionTwo(int index) const { return lengthFractionTwo[index]; }
  void SelectedKiWaTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout);

 private:
  int numberElements;
  int * timeSeriesElements;
  double * lengthFractionOne;
  double * lengthFractionTwo;
};

SelectedKiWaTimeSeriesElements::SelectedKiWaTimeSeriesElements():
  numberElements(0)
{
}
     
SelectedKiWaTimeSeriesElements::~SelectedKiWaTimeSeriesElements()
{     
}

void SelectedKiWaTimeSeriesElements::SelectedKiWaTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout)
{
  char fileName[80], buffer[256];
  int i, numElements, elementId;
  double lthFracOne, lthFracTwo;

  /*  cout << " File with landscape elements selected for KiWa time series output: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream fileTimeSeries(fileName);  // Open for reading
  if (fileTimeSeries == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  fileTimeSeries.ignore(100,':');
  fileTimeSeries >> numElements;
  SetNumberElements(numElements);
  timeSeriesElements = new int [numElements];
  lengthFractionOne = new double [numElements];
  lengthFractionTwo = new double [numElements];
  fileTimeSeries.ignore(256,'\n');
  fileTimeSeries.getline(buffer, 256);
  for (i=0; i<numElements; i++) {
    fileTimeSeries >> elementId >> lthFracOne >> lthFracTwo;
    if (lthFracOne<0.0 || lthFracOne>1.0) lthFracOne=0.5;
    if (lthFracTwo<0.0 || lthFracTwo>1.0) lthFracTwo=1.0;
    SetKiWaTimeSeriesElement(i, elementId);
    SetLengthFractionOne(i, lthFracOne);
    SetLengthFractionTwo(i, lthFracTwo);
    fileTimeSeries.ignore(256,'\n');
  }
  fout << endl << "Number of landscape elements selected for KiWa time series output: " << GetNumberElements() << endl;
  for (i=0; i<GetNumberElements(); i++) {
    fout << GetKiWaTimeSeriesElement(i) << "  " << GetLengthFractionOne(i) << "  " << GetLengthFractionTwo(i) << endl;
  }
  fout << endl;
}


//class SelectedHbvTimeSeriesElements
class SelectedHbvTimeSeriesElements
{
 public:
  SelectedHbvTimeSeriesElements();
  ~SelectedHbvTimeSeriesElements();
  void SetNumberElements(int value) { numberElements = value; }
  int GetNumberElements() const { return numberElements; }
  void SetHbvTimeSeriesElement(int index, int value) { timeSeriesElements[index] = value; }
  int GetHbvTimeSeriesElement(int index) const { return timeSeriesElements[index]; }
  void SetEffPor(int index, double value) { effectivePorosity[index] = value; }
  double GetEffPor(int index) const { return effectivePorosity[index]; }
  void SetGroundWaterRef(int index, double value) { groundWaterReference[index] = value; }
  double GetGroundWaterRef(int index) const { return groundWaterReference[index]; }
  void SelectedHbvTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout);

 private:
  int numberElements;
  int * timeSeriesElements;
  double * effectivePorosity;
  double * groundWaterReference;
};

SelectedHbvTimeSeriesElements::SelectedHbvTimeSeriesElements():
  numberElements(0)
{
}
     
SelectedHbvTimeSeriesElements::~SelectedHbvTimeSeriesElements()
{     
}

void SelectedHbvTimeSeriesElements::SelectedHbvTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout)
{
  char fileName[80], buffer[256];
  int i, numElements, elementId;
  double effPor, gwRef;

  /*  cout << " File with landscape elements selected for HBV time series output: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream fileTimeSeries(fileName);  // Open for reading
  if (fileTimeSeries == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  fileTimeSeries.ignore(100,':');
  fileTimeSeries >> numElements;
  SetNumberElements(numElements);
  timeSeriesElements = new int [numElements];
  effectivePorosity = new double [numElements];
  groundWaterReference = new double [numElements];
  fileTimeSeries.ignore(256,'\n');
  fileTimeSeries.getline(buffer, 256);
  for (i=0; i<numElements; i++) {
    fileTimeSeries >> elementId >> gwRef >> effPor;
    SetHbvTimeSeriesElement(i, elementId);
    SetGroundWaterRef(i, gwRef);
    SetEffPor(i, effPor);
    fileTimeSeries.ignore(256,'\n');
  }
  fout << "Number of landscape elements selected for HBV time series output: " << GetNumberElements() << endl;
  for (i=0; i<GetNumberElements(); i++) {
    fout << GetHbvTimeSeriesElement(i) << "  " << GetGroundWaterRef(i) << "  " << GetEffPor(i) << endl;
  }
  fout << endl;
}


//class SelectedWaterCourseTimeSeriesElements
class SelectedWaterCourseTimeSeriesElements
{
 public:
  SelectedWaterCourseTimeSeriesElements();
  ~SelectedWaterCourseTimeSeriesElements();
  void SetNumberElements(int value) { numberElements = value; }
  int GetNumberElements() const { return numberElements; }
  void SetWaterCourseTimeSeriesElement(int index, int value) { timeSeriesElements[index] = value; }
  int GetWaterCourseTimeSeriesElement(int index) const { return timeSeriesElements[index]; }
  void SelectedWaterCourseTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout);

 private:
  int numberElements;
  int * timeSeriesElements;
};

SelectedWaterCourseTimeSeriesElements::SelectedWaterCourseTimeSeriesElements():
  numberElements(0)
{
}
     
SelectedWaterCourseTimeSeriesElements::~SelectedWaterCourseTimeSeriesElements()
{     
}

void SelectedWaterCourseTimeSeriesElements::SelectedWaterCourseTimeSeriesElementsInput(ifstream &fileControl, ofstream &fout)
{
  char fileName[80], buffer[256];
  char ch;
  int i, j, numElements, elementId;

  /*  cout << " File with watercourse/catchment elements selected for time series output: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream fileTimeSeries(fileName);  // Open for reading
  if (fileTimeSeries == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  fileTimeSeries.ignore(100,':');
  fileTimeSeries >> numElements;
  SetNumberElements(numElements);
  timeSeriesElements = new int [numElements];
  fileTimeSeries.ignore(256,'\n');
  for (i=0; i<numElements; i++) {
    fileTimeSeries >> j >> ch >> elementId;
    if (j != i) {
      cout << endl << "Error reading file " << fileName << "\t" << i << "\t" << j << endl;
      exit (1);
    }
    SetWaterCourseTimeSeriesElement(i, elementId);
    fileTimeSeries.ignore(256,'\n');
  }
  fout << "Number of watercourse/catchment elements selected for time series output: " << GetNumberElements() << endl;
  for (i=0; i<GetNumberElements(); i++) {
    fout << GetWaterCourseTimeSeriesElement(i) << endl;
  }
  fout << endl;
}


// class TotalReservoirStorage
class TotalReservoirStorage
{
 public:
  TotalReservoirStorage();
  ~TotalReservoirStorage();
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  void AllocateTotalReservoirStorage(int numberTimeSteps); 
  void SetInitialTotalReservoirStorage();
  double GetInitialTotalReservoirStorage() const { return initialTotalReservoirStorage; };
  void SetTotalReservoirStorage(int index, double value) { totalReservoirStorage[index] = value; }
  double GetTotalReservoirStorage(int index) const { return totalReservoirStorage[index]; }
  
 private:
  ParametersCommon *commonPar;
  double initialTotalReservoirStorage;
  double *totalReservoirStorage;
};

TotalReservoirStorage::TotalReservoirStorage():
  initialTotalReservoirStorage(0)
{
  SetCommonPar(0);
}
     
TotalReservoirStorage::~TotalReservoirStorage()
{     
}

void TotalReservoirStorage::SetInitialTotalReservoirStorage()
{
  initialTotalReservoirStorage = commonPar->GetINITIAL_TOTAL_RESERVOIR();
}

void TotalReservoirStorage::AllocateTotalReservoirStorage(int numberTimeSteps) 
{
  int i;
  totalReservoirStorage = new double [numberTimeSteps];
  for (i=0; i<numberTimeSteps; i++) totalReservoirStorage[i]=missingData;
}


// class DateTimeInformation
class DateTimeInformation
{
 public:
  DateTimeInformation();
  ~DateTimeInformation();
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  void SetInitialTimeSteps(int value) { initialTimeSteps = value; }
  int GetInitialTimeSteps() { return initialTimeSteps; }
  void SetNumberTimeSteps(int value) { numberTimeSteps = value; }
  int GetNumberTimeSteps() { return numberTimeSteps; }
  void SetStartModelTime(DateTime datetime) { startModelTime = datetime; }
  DateTime GetStartModelTime() { return startModelTime; }
  void SetStartSimulationTime(DateTime datetime) { startSimulationTime = datetime; }
  DateTime GetStartSimulationTime() { return startSimulationTime; }
  void SetEndSimulationTime(DateTime datetime) { endSimulationTime = datetime; }
  DateTime GetEndSimulationTime() { return endSimulationTime; }
  void ReadDateTimeInformation(ifstream &fileControl);
  void CalculateTimeStepInformation();
  void WriteDateTimeInformation(ofstream &fout);

 private:
  ParametersCommon *commonPar;
  int initialTimeSteps;
  int numberTimeSteps;
  DateTime startModelTime;
  DateTime startSimulationTime;
  DateTime endSimulationTime;
};

DateTimeInformation::DateTimeInformation():
  numberTimeSteps(0),
  initialTimeSteps(0)
{
  SetCommonPar(0);
}

DateTimeInformation::~DateTimeInformation()
{
}

void DateTimeInformation::ReadDateTimeInformation(ifstream &fileControl)
{
  int day, mth, year, hour, minute;
  /*  cout << " Start model date and time (day, month, year, hour, minute): ";
      cin >> day >> mth >> year >> hour >> minute;*/
  fileControl.ignore(100,':');
  fileControl >> day >> mth >> year >> hour >> minute;
  fileControl.ignore(256,'\n');
  DateTime startModel(year,mth,day,hour,minute,0);
  if (startModel.legal() != 1) {
    cout << endl << " Not legal date start model and time " << endl << endl;
    exit(1);
  }
  /*  cout << endl;
      cout << " Start simulation date and time (day, month, year, hour, minute): ";
      cin >> day >> mth >> year >> hour >> minute;*/
  fileControl.ignore(100,':');
  fileControl >> day >> mth >> year >> hour >> minute;
  fileControl.ignore(256,'\n');
  DateTime startSimulation(year,mth,day,hour,minute,0);
  if (startSimulation.legal() != 1) {
    cout << endl << " Not legal start simulation date and time " << endl << endl;
    exit(1);
  }
  /*  cout << endl;
      cout << " End simulation date and time (day, month, year, hour, minute): ";
      cin >> day >> mth >> year >> hour >> minute;*/
  fileControl.ignore(100,':');
  fileControl >> day >> mth >> year >> hour >> minute;
  fileControl.ignore(256,'\n');
  DateTime endSimulation(year,mth,day,hour,minute,0);
  if (endSimulation.legal() != 1) {
    cout << endl << " Not legal end simulation date and time " << endl << endl;
    exit(1);
  }
  cout << endl;
  SetStartModelTime(startModel);
  SetStartSimulationTime(startSimulation);
  SetEndSimulationTime(endSimulation);
  if (startModelTime > startSimulationTime || startSimulationTime > endSimulationTime) {
    cout << endl << " startModelTime " << startModelTime << "    startSimulationTime " << startSimulationTime;
    cout << "    endSimulationTime " << endSimulationTime << endl << endl;
    exit(1);
  }
}

void DateTimeInformation::CalculateTimeStepInformation()
{
  int initTime, numTime;
  initTime = (int)((startSimulationTime-startModelTime)/(double)commonPar->GetSECONDS_TIMESTEP());
  numTime = 1+(int)((endSimulationTime-startSimulationTime)/(double)commonPar->GetSECONDS_TIMESTEP());
  //  initialTimeSteps = (int)(startSimulationTime.date2jday()-startModelTime.date2jday());
  //  numberTimeSteps = (int)(endSimulationTime.date2jday()-startSimulationTime.date2jday()+1);
  SetInitialTimeSteps(initTime);
  SetNumberTimeSteps(numTime);
}

void DateTimeInformation::WriteDateTimeInformation(ofstream &fout)
{
  fout << "Date and time information: \n";
  fout << "Start model time      " << GetStartModelTime() << endl; 
  fout << "Start simulation time " << GetStartSimulationTime() << endl; 
  fout << "End simulation time   " << GetEndSimulationTime() << endl; 
  fout << "Initial time steps    " << GetInitialTimeSteps() << endl; 
  fout << "Number time steps     " << GetNumberTimeSteps() << endl; 
  fout << endl;
}


// class MeteorologicalStations
class MeteorologicalStations
{
 public:
  MeteorologicalStations();
  ~MeteorologicalStations();
  void SetNumPrecStations(int value) { numberPrecStations = value; }
  int GetNumPrecStations() { return numberPrecStations; }
  void SetNumTempStations(int value) { numberTempStations = value; }
  int GetNumTempStations() { return numberTempStations; }
  int GetStationNumber(int k)  { return stationNumber[k]; }
  double GetStationCoordX(int k)  { return stationCoordX[k]; }
  double GetStationCoordY(int k)  { return stationCoordY[k]; }
  double GetStationAltitude(int k)  { return stationAltitude[k]; }
  void SetMeteorologicalStations(ifstream &fileControl, ofstream &fout);

 private:
  int numberPrecStations;
  int numberTempStations;
  int * stationNumber;
  double * stationCoordX;
  double * stationCoordY;
  double * stationAltitude;
};

MeteorologicalStations::MeteorologicalStations():
  numberPrecStations(0),
  numberTempStations(0)
{
}

MeteorologicalStations::~MeteorologicalStations()
{
}

void MeteorologicalStations::SetMeteorologicalStations(ifstream &fileControl, ofstream &fout) 
{
  int i, numberPrecStations, numberTempStations;
  char stationType;
  char fileName[80];
  /*  cout << " File with meteorological stations: ";
      cin >> fileName;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  //  cout << fileName << endl;
  ifstream finMetSta(fileName);  // Open for reading
  if (finMetSta == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finMetSta.ignore(100,':'); finMetSta >> numberPrecStations;
  finMetSta.ignore(100,':'); finMetSta >> numberTempStations;
  SetNumPrecStations(numberPrecStations);
  SetNumTempStations(numberTempStations);
  stationNumber = new int [numberPrecStations+numberTempStations];
  stationCoordX = new double [numberPrecStations+numberTempStations];
  stationCoordY = new double [numberPrecStations+numberTempStations];
  stationAltitude = new double [numberPrecStations+numberTempStations];
  for (i=0; i<numberPrecStations; i++) {
    finMetSta >> stationType >> stationNumber[i] >> stationCoordX[i] >> stationCoordY[i] >> stationAltitude[i];
    finMetSta.ignore(256,'\n');
    if (stationType != 'P' && stationType != 'p') {
      cout << endl << " Not legal station type : " << stationType << " , Station number = " << i << endl << endl;
      exit(1);
    }
  }
  for (i=0; i<numberTempStations; i++) {
    finMetSta >> stationType >> stationNumber[numberPrecStations+i] >> stationCoordX[numberPrecStations+i] 
        >> stationCoordY[numberPrecStations+i] >> stationAltitude[numberPrecStations+i];
    finMetSta.ignore(256,'\n');
    if (stationType != 'T' && stationType != 't') {
      cout << endl << " Not legal station type : " << stationType  
           << " station number = " << numberPrecStations+i << endl << endl;
      exit(1);
    }
  }
  finMetSta.close();
  fout << endl << "Meteorological stations: \n";
  fout << GetNumPrecStations() << endl;
  fout << GetNumTempStations() << endl;
  fout << "Prec.    " << endl; 
  for (i=0; i<GetNumPrecStations(); i++) {
    fout << GetStationNumber(i) << "  ";
    fout << GetStationCoordX(i) << "  ";
    fout << GetStationCoordY(i) << "  ";
    fout << GetStationAltitude(i) << endl;
  }
  fout << "Temp.    " << endl;
  for (i=0; i<GetNumTempStations(); i++) {
    fout << GetStationNumber(GetNumPrecStations()+i) << "  ";
    fout << GetStationCoordX(GetNumPrecStations()+i) << "  ";
    fout << GetStationCoordY(GetNumPrecStations()+i) << "  ";
    fout << GetStationAltitude(GetNumPrecStations()+i) << endl;
  }
  fout << endl;
}


// class InputTimeSeries
class InputTimeSeries
{
 public:
  InputTimeSeries(int numberRows, int numberColums, DateTime firstTime, DateTime lastTime, int secondsPerTimeStep);
  ~InputTimeSeries();
  int GetNumberTimeSteps() const { return timeSteps; }
  int GetNumberInputSeries() const { return numberSeries; }
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  DateTime GetDateTime(int i) const { return datetime[i]; }
  void SetInput(ifstream &fin);
  double GetInput(int i, int j) const { return inputArray[i*numberSeries+j]; }
  void WriteInput();

 private:
  ParametersCommon *commonPar;
  int timeSteps;
  int numberSeries;
  DateTime * datetime;
  double * inputArray;
};

InputTimeSeries::InputTimeSeries(int numberRows, int numberColums, DateTime firstTime, DateTime lastTime, int secondsPerTimeStep):
  timeSteps(numberRows),
  numberSeries(numberColums)
{
  int i;
  SetCommonPar(0);
  datetime = new DateTime [timeSteps];
  inputArray = new double [timeSteps*numberSeries];
  for (i=0; i<timeSteps; i++) {
    datetime[i] = firstTime + i*secondsPerTimeStep;
  }
  for (i=0; i<timeSteps*numberSeries; i++) {
    inputArray[i] = missingData;
    //    cout << i << " " << inputArray[i] << "  ";
  }
  if (datetime[0] != firstTime || datetime[timeSteps-1] != lastTime) {
    cout << endl << " DateTime error during initialisation of InputTimeSeries array: " <<
      datetime[0] << "  " << firstTime  << "  " << datetime[timeSteps-1]  << "  " << lastTime << endl << endl;
    exit(1);
  }
  //  cout << datetime[0] << "  " << firstTime  << "  " << datetime[timeSteps-1]  << "  " << lastTime << endl << endl;
}

InputTimeSeries::~InputTimeSeries()
{
}

void InputTimeSeries::SetInput(ifstream &fin)
{
  char ch;
  char buffer[256];
  int i, j;
  int date, time, year, mth, day, hour, min;
  double * tempArray = new double [numberSeries];
  fin.getline(buffer, 256);
  fin.getline(buffer, 256);
  i = 0;
  while (fin >> date >> ch >> time) {
    for (j=0;j<numberSeries;j++) {
      fin >> tempArray[j];
    }
    year = date/10000;
    mth = (date%10000)/100;
    day = date%100;
    hour = time/100;
    min = time%100;
    //    cout << year << " " << mth << " " << day << " " << hour << " " << min << " " << endl;
    while ((datetime[i].getYear()<year ||
       (datetime[i].getYear()==year && datetime[i].getMonth()<mth) ||
       (datetime[i].getYear()==year && datetime[i].getMonth()==mth && datetime[i].getDay()<day) ||
       (datetime[i].getYear()==year && datetime[i].getMonth()==mth && datetime[i].getDay()==day && datetime[i].getHour()<hour) ||
       (datetime[i].getYear()==year && datetime[i].getMonth()==mth && datetime[i].getDay()==day && 
        datetime[i].getHour()==hour && datetime[i].getMinute()<min)) &&
        i<timeSteps-1) {
      i++;
      //      cout << " i = " << i << endl;
    }
    //    cout << datetime[i].getYear() << " " << datetime[i].getMonth() << " " << datetime[i].getDay() << " " << 
    //      datetime[i].getHour() << " " << datetime[i].getMinute() << endl;
    if (datetime[i].getYear()==year && datetime[i].getMonth()==mth && datetime[i].getDay()==day && 
        datetime[i].getHour()==hour && datetime[i].getMinute()==min) {
      //      cout << " data found  " << i << " " << year << " " << mth << " " << day << " " << hour << " " << min << " " << endl;
      for (j=0;j<numberSeries;j++) {
        inputArray[i*numberSeries+j] = tempArray[j];
      }
    }
  }
  delete [] tempArray;
}

void InputTimeSeries::WriteInput()
{
  FILE *fp_out;
  int i, j;
  if ((fp_out = fopen("input_out.txt", "w")) == NULL ) {
    printf("\n File input_out.txt not found!\n\n");
    exit(1);
  }
  for (i=0; i<timeSteps; i++) {
    fprintf(fp_out,"%04d%02d%02d/%02d%02d",GetDateTime(i).getYear(),GetDateTime(i).getMonth(),GetDateTime(i).getDay(),
            GetDateTime(i).getHour(),GetDateTime(i).getMinute());
    for (j=0;j<numberSeries;j++){
      fprintf(fp_out,"%15.5f",GetInput(i,j));
    }
    fprintf(fp_out,"\n");
  }
  fclose(fp_out);
}


// class InputElement
class InputElement
{
 public:
  InputElement(int value);
  ~InputElement();
  int GetNumberValues() const { return numberValues; }
  void SetInput(int i, double value) { inputArray[i] = value; }
  double GetInput(int i) const { return inputArray[i]; }

 private:
  int numberValues;
  double * inputArray;
};

InputElement::InputElement(int value):
  numberValues(value)
{
  int i;
  inputArray = new double [numberValues];
  for (i=0; i<numberValues; i++) inputArray[i] = missingData;
}

InputElement::~InputElement()
{
}


// class LakeWaterBalance
class LakeWaterBalance
{
 public:
  LakeWaterBalance();
  ~LakeWaterBalance();
  void WaterBalance(int timeStep, double upLandAccumulatedDischarge);
  void SetLakeValues(double temperature, double level);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetLakeEvap() const;
  double GetRunoff() const;
  double GetWaterLevelChange() const;
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }

 private:
  ParametersCommon *commonPar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedElement *landScapeElement;
  double precipitation;                       /*  Precipitation (m/timestep)  */
  double temp;                                /*  Air temperature (deg. C)   */
  double lakeTemp;                            /*  Lake temperature (deg. C)  */
  double lakeEvaporation;                     /*  Evaporation (m/timestep)  */
  double waterLevel;                          /*  Lake water level (m)  */
  double waterLevelChange;                    /*  Lake water level change (m)  */
  double runoff;                              /*  Runoff (m/timestep)  */
  double discharge;                           /*  Lake outflow (m3/s)  */
};

LakeWaterBalance::LakeWaterBalance():
  precipitation(0.0),
  temp(0.0),
  waterLevelChange(0.0),
  lakeEvaporation(0.0),
  runoff(0.0),
  discharge(0.0)
{
  SetCommonPar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

LakeWaterBalance::~LakeWaterBalance()
{
}

void LakeWaterBalance::SetLakeValues(double temperature, double level)
{
  lakeTemp = temperature;
  waterLevel = level;
}

double LakeWaterBalance::GetPrecipitation() const
{
  return precipitation;
}

double LakeWaterBalance::GetTemperature() const
{
  return temp;
}

double LakeWaterBalance::GetLakeEvap() const
{
  return lakeEvaporation;
}

double LakeWaterBalance::GetRunoff() const
{
  return runoff;
}

double LakeWaterBalance::GetWaterLevelChange() const
{
  return waterLevelChange;
}


// class Vegetation
class Vegetation
{
 public:
  Vegetation();
  ~Vegetation();
  void WaterBalance(int timeStep);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetInterceptionLoss() const;
  double GetThroughFall() const;
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }
  double GetDryPeriod() { return dryPeriod; }

 private:
  ParametersCommon *commonPar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedElement *landScapeElement;
  double precipitation;                    /*  Precipitation (m/timestep)  */
  double temp;                             /*  Air temperature (deg. C)  */
  double potev;                            /*  Potential evapotranspiration (m/timestep)  */
  double prevInterception;                 /*  Interception from previous time step (m)  */
  double interceptionStore;                /*  Interception store (m)  */
  double interceptionLoss;                 /*  Interception loss from vegetation (m)  */
  double throughFall;                      /*  Throughfall (m/timestep)  */
  double wetPeriod, dryPeriod;             /*  Length of wet and dry period during evapotranspiration (fraction of timestep)  */
};

Vegetation::Vegetation():
  precipitation(0.0),
  temp(0.0),
  potev(0.0),
  prevInterception(0.0),
  interceptionStore(0.0),
  interceptionLoss(0.0),
  throughFall(0.0),
  wetPeriod(0.0),
  dryPeriod(0.0)
{
  SetCommonPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

Vegetation::~Vegetation()
{
}

double Vegetation::GetPrecipitation() const
{
  return precipitation;
}

double Vegetation::GetTemperature() const
{
  return temp;
}

double Vegetation::GetInterceptionLoss() const
{
  return interceptionLoss;
}

double Vegetation::GetThroughFall() const
{
  return throughFall;
}


// class Snow
class Snow
{
 public:
  Snow();
  ~Snow();
  void WaterBalance(int timeStep, double waterInput);
  void SetSnowStore(double value);
  double GetSnowCoverFraction() const;
  double GetSnowStore() const;
  double GetSnowWaterEquivalentChange() const;
  double GetMeltWater() const;
  double GetWaterOutput() const;
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }

 private:
  ParametersCommon *commonPar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedElement *landScapeElement;
  double temp;                               /*  Air temperature (deg. C)  */
  double snowStore;                          /*  Snow store (m)  */
  double meltWater;                          /*  Meltwater in snow (m)  */
  double waterOutput;                        /*  Output of meltwater from snow store (m/timestep)  */
  double snowWaterEquivalentChange;          /*  Change of snow water equivalent (m)  */
  double snowCoverFraction;                  /*  Fraction of area covered by snow  */
  double distSnowStore[numberSnowClasses];   /*  Distributed snow store (m)  */
  double distMeltWater[numberSnowClasses];   /*  Distributed meltwater in snow (m)  */
  double distWaterOutput[numberSnowClasses]; /*  Output of meltwater from snow store (m/timestep)  */
};

Snow::Snow():
  temp(0.0),
  meltWater(0.0),
  snowWaterEquivalentChange(0.0),
  waterOutput(0.0),
  snowCoverFraction(0.0)
{
  int i;
  for (i=0; i<numberSnowClasses; i++) {
    distSnowStore[i]=0.0;
    distMeltWater[i]=0.0;
    distWaterOutput[i]=0.0;
  }
  SetCommonPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

Snow::~Snow()
{
}

void Snow::SetSnowStore(double value)
{
  snowStore = value;
  //  meltWater = value*landSurfacePar->GetMAX_REL();
  for (int i=0; i<numberSnowClasses; i++) {
    distSnowStore[i] = value;
    //    distMeltWater[i] = value*landSurfacePar->GetMAX_REL();
  }
}

double Snow::GetSnowCoverFraction() const
{ 
return snowCoverFraction; 
}

double Snow::GetSnowStore() const
{
  return snowStore;
}

double Snow::GetMeltWater() const
{
  return meltWater;
}

double Snow::GetSnowWaterEquivalentChange() const
{
  return snowWaterEquivalentChange;
}

double Snow::GetWaterOutput() const
{
  return waterOutput;
}


// class GlacierSurface
class GlacierSurface
{
 public:
  GlacierSurface();
  ~GlacierSurface();
  void WaterBalance(int timeStep, int initialTimeSteps);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetIceMelt() const;
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }

 private:
  ParametersCommon *commonPar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedElement *landScapeElement;
  double precipitation;                    /*  Precipitation (m/timestep)  */
  double temp;                             /*  Air temperature (deg. C)  */
  double iceMelt;                          /*  Meltwater from ice (m/timestep)  */
};

GlacierSurface::GlacierSurface():
  precipitation(0.0),
  temp(0.0),
  iceMelt(0.0)
{
  SetCommonPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

GlacierSurface::~GlacierSurface()
{
}

double GlacierSurface::GetPrecipitation() const
{
  return precipitation;
}

double GlacierSurface::GetTemperature() const
{
  return temp;
}

double GlacierSurface::GetIceMelt() const
{
  return iceMelt;
}


// class GlacierIce                        /*  Currently not in use !  */
class GlacierIce
{
 public:
  GlacierIce();
  ~GlacierIce();
  void WaterBalance(int timeStep, double waterInput);
  double GetWaterOutput() const;
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }

 private:
  ParametersCommon *commonPar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedElement *landScapeElement;
  double waterStored;
  double waterOutput;
};

GlacierIce::GlacierIce():
  waterStored(0.0),
  waterOutput(0.0)
{
  SetCommonPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

GlacierIce::~GlacierIce()
{
}

void GlacierIce::WaterBalance(int timeStep, double waterInput)
{
  waterStored=waterInput;
  waterOutput=waterStored;
  waterStored=waterStored-waterOutput;
}

double GlacierIce::GetWaterOutput() const
{
  return waterOutput;
}


// class HBV
class HBV
{
 public:
  HBV();
  ~HBV();
  void SetInitialHbvValues();
  void SetSubSurfaceHbvStore(double sm, double uz, double lz);
  // Algorithm to be performed in case: no input to landscape element from upstream elements
  //  void WaterBalance(int timeStep, double waterInput, double dryPeriod);
  // Algorithm to be performed in case: input to landscape element from upstream elements
  void WaterBalance(int timeStep, double waterInput, double snowCoverFraction, double dryPeriod, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge);
  double GetSoilMoisture() const;
  double GetSoilMoistureDeficit() const;
  double GetPercSoilUpper() const;
  double GetUpperZone() const;
  double GetLowerZone() const;
  double GetTranspSoilEvap() const;
  double GetLowerRunoff() const;
  double GetUpperRunoff() const;
  double GetRunoff() const;
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  void SetSubSurfaceHbvPar(ParametersSubSurfaceHbv *parObj) { subSurfacePar = parObj; }
  ParametersSubSurfaceHbv *GetSubSurfaceHbvPar() const { return subSurfacePar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }

 private:
  ParametersCommon *commonPar;
  ParametersSubSurfaceHbv *subSurfacePar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedElement *landScapeElement;
  double temp;                                /*  Temperature (deg. C)  */
  double soilMoisture;                        /*  Soil moisture content (m)  */ 
  double percSoilUpper;                       /*  Percolation from soil moisture zone to upper zone (m/timestep)  */ 
  double upperZone;                           /*  Upper groundwater zone water content (m)  */
  double lowerZone;                           /*  Lower groundwater zone water content (m)  */
  double lowerRunoff;                         /*  Runoff from lower layer (m/timestep)  */
  double upperRunoff;                         /*  Runoff from upper layer (m/timestep)  */
  double transpSoilEvap;                      /*  Water lost from subsurface by evapotranspiration (m)  */
  double runoff;                              /*  Runoff (m/timestep)  */
};

HBV::HBV():
  temp(0.0),
  soilMoisture(0.0),
  percSoilUpper(0.0),
  upperZone(0.0),
  lowerZone(0.0),
  transpSoilEvap(0.0),
  lowerRunoff(0.0),
  upperRunoff(0.0),
  runoff(0.0)
{
  SetCommonPar(0);
  SetSubSurfaceHbvPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
}

HBV::~HBV()
{
}

void HBV::SetInitialHbvValues()
{
  soilMoisture = commonPar->GetINITIAL_SOIL_MOISTURE();
  upperZone = commonPar->GetINITIAL_UPPER_ZONE();
  lowerZone = commonPar->GetINITIAL_LOWER_ZONE();
}

void HBV::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
  soilMoisture = sm;
  upperZone = uz;
  lowerZone = lz;
}

double HBV::GetSoilMoisture() const
{
  return soilMoisture;
}

double HBV::GetSoilMoistureDeficit() const
{
  return subSurfacePar->GetFC()-soilMoisture;
}

double HBV::GetPercSoilUpper() const
{
  return percSoilUpper;
}

double HBV::GetUpperZone() const
{
  return upperZone;
}

double HBV::GetLowerZone() const
{
  return lowerZone;
}

double HBV::GetTranspSoilEvap() const
{
  return transpSoilEvap;
}

double HBV::GetRunoff() const
{
  return runoff;
}

double HBV::GetLowerRunoff() const
{
  return lowerRunoff;
}

double HBV::GetUpperRunoff() const
{
  return upperRunoff;
}


// class KinematicWave
class KinematicWave
{
 public:
  KinematicWave();
  ~KinematicWave();
  void SetInitialKwaValues();
  void WaterBalance(int timeStep, double waterInput, double snowCoverFraction, double dry_period, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge);
  void kinematic_wave_with_lateral_inflow (double, double *, double *, double *, double *, double *, double *, double *, int *,
                                           double, double, double, double, double, double, double, double);
  void kinematic_wave_without_lateral_inflow (double *, double *, double *, double *, double *, double *, double *, int *, 
                                              double, double, double, double, double, double, double, double);
  void KiWaGroundWaterTable(int elementId, int timeStep);
  void KiWaSoilMoisture(int elementId, int timeStep);
  double GetSoilMoisture(double lengthFraction) const;
  double GetGroundWaterDepth(double lengthFraction) const;
  double GetTranspSoilEvap() const;
  double GetLowerRunoff() const;
  double GetUpperRunoff() const;
  double GetRunoff() const;
  void SetSelectedKiWaHillslopeElements(SelectedKiWaHillslopeElements *object) { selectedKiWaElements = object; }
  SelectedKiWaHillslopeElements *GetSelectedKiWaHillslopeElements() const { return selectedKiWaElements; }
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  void SetLandSurfacePar(ParametersLandSurface *parObj) { landSurfacePar = parObj; }
  ParametersLandSurface *GetLandSurfacePar() { return landSurfacePar; }
  void SetKiWaPar(ParametersKiWa *parObj) { kiWaPar = parObj; }
  ParametersKiWa *GetKiWaPar() const { return kiWaPar; }
  //  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
  //  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
  void SetInputElement(InputElement *inElementObj) { inElement = inElementObj; }
  InputElement *GetInputElement() const { return inElement; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }
  void SetDateTimeInfo(DateTimeInformation *dateTimeObj) { dateTimeInfo = dateTimeObj; }
  DateTimeInformation *GetDateTimeInfo() const { return dateTimeInfo; }

 private:
  SelectedKiWaHillslopeElements *selectedKiWaElements;
  ParametersCommon *commonPar;
  ParametersKiWa *kiWaPar;
  ParametersLandSurface *landSurfacePar;
  //  InputTimeSeries *inTimeSeries;
  InputElement *inElement;
  DistributedElement *landScapeElement;
  DateTimeInformation *dateTimeInfo;
  int u_high;                                 /*  Number of characteristic curves in upper layer  */
  double temp;                                /*  Temperature (deg. C)  */
  double vap_def;                             /*  Vapor pressure deficit (hPa)  */
  double mean_perc;                           /*  Mean percolation along hillslope (m/time_step)  */
  double field_capacity;                      /*  Field capacity  */
  double teta;                                /*  Volumetric water content */
  double gw_h;                                /*  Depth to groundwater table (m)  */
  double evaporated_volume;                   /*  Volume of water lost from lower layer by evapotranspiration (m)  */
  double evol_upper;                          /*  Volume of water lost from upper layer by evapotranspiration (m)  */
  double volume_root;                         /*  Volume of water removed from root zone by evapotranspiration (m)  */
  double def_par;                             /*  Fraction of actual evapotranspiration removed from soil moisture  */
  double smdef[numberCharacteristic+1];       /*  Soil moisture deficit (m) ( <= 0 )  */ 
  double perc[numberCharacteristic+1];        /*  Volume of water percolating to saturated zone (m)  */
  double len_coord[numberCharacteristic+1];   /*  Length-coordinate along characteristic curves in lower layer (m)  */
  double fixed_length[numberCharacteristic+1];/*  Fixed length-coordinate along hillslope (m)  */
  double sat_depth[numberCharacteristic+1];   /*  Depth of saturated zone along characteristic curves in lower layer (m)  */
  double fixed_sat[numberCharacteristic+1];   /*  Saturated depth at fixed length coordinates  */
  double upp_tim[numberCharacteristic+1];     /*  Initial time within time step of characteristic curves in upper layer  */
  double upp_dep[numberCharacteristic+1];     /*  Depth along characteristic curves in upper layer  */
  double upp_len[numberCharacteristic+1];     /*  Length coordinate along characteristic curves in upper layer  */
  double actev;                               /*  Actual evapotranspiration (m/time_step)  */
  double actev_loss;                          /*  Actual evapotranspiration loss from soil (m)  */
  double transpSoilEvap;                      /*  Water lost from subsurface by evapotranspiration (m)  */
  double lower_runoff;                        /*  Runoff from lower layer (mm/time_step)  */
  double upper_runoff;                        /*  Runoff from upper layer (mm/time_step)  */
  double runoff;                              /*  Runoff (m/timestep)  */
};

KinematicWave::KinematicWave():
  temp(0.0),
  vap_def(0.0),
  u_high(0),
  mean_perc(0.0),
  field_capacity(0.0),
  teta(0.0),
  gw_h(0.0),
  evaporated_volume(0.0),
  evol_upper(0.0),
  volume_root(0.0),
  def_par(0.0),
  actev(0.0),
  actev_loss(0.0),
  transpSoilEvap(0.0),
  lower_runoff(0.0),
  upper_runoff(0.0),
  runoff(0.0)
{
  SetSelectedKiWaHillslopeElements(0);
  SetCommonPar(0);
  SetKiWaPar(0);
  SetLandSurfacePar(0);
  //  SetInputTimeSeries(0);
  SetInputElement(0);
  SetLandScapeElement(0);
  SetDateTimeInfo(0);
}

KinematicWave::~KinematicWave()
{
}

void KinematicWave::SetInitialKwaValues()
{
}

double KinematicWave::GetTranspSoilEvap() const
{
  return transpSoilEvap;
}

double KinematicWave::GetRunoff() const
{
  return runoff;
}

double KinematicWave::GetLowerRunoff() const
{
  return lower_runoff;
}

double KinematicWave::GetUpperRunoff() const
{
  return upper_runoff;
}


// class Lake
class Lake
{
 public:
  Lake();
  ~Lake();
  void SetLakeWaterBalance(LakeWaterBalance *theLakeWaterBalance) { ptrLakeWaterBalance = theLakeWaterBalance; }
  LakeWaterBalance *GetLakeWaterBalance() const { return ptrLakeWaterBalance; }
  void SetAreaFraction(double value) { areaFraction = value; }
  double GetAreaFraction() const { return areaFraction; }
  void WaterBalance(int timeStep, double upLandAccumulatedDischarge);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetLakeEvap() const;
  double GetRunoff() const;
  double GetLakeStorageChange() const;

 private:
  LakeWaterBalance *ptrLakeWaterBalance;
  double areaFraction;
};

Lake::Lake()
{
  SetAreaFraction(0.0);
}

Lake::~Lake()
{
}

void Lake::WaterBalance(int timeStep, double upLandAccumulatedDischarge)
{
  GetLakeWaterBalance()->WaterBalance(timeStep, upLandAccumulatedDischarge);
}

double Lake::GetPrecipitation() const
{
  double precipitation=GetLakeWaterBalance()->GetPrecipitation();
  return precipitation*GetAreaFraction()/100.0;
}

double Lake::GetTemperature() const
{
  double temperature=GetLakeWaterBalance()->GetTemperature();
  return temperature*GetAreaFraction()/100.0;
}

double Lake::GetLakeEvap() const
{
  double lakeEvap=GetLakeWaterBalance()->GetLakeEvap();
  return lakeEvap*GetAreaFraction()/100.0;
}

double Lake::GetRunoff() const
{
  double runoff=GetLakeWaterBalance()->GetRunoff();
  return runoff*GetAreaFraction()/100.0;
}

double Lake::GetLakeStorageChange() const
{
  double storageChange=GetLakeWaterBalance()->GetWaterLevelChange()*GetAreaFraction()/100.0;
  return storageChange;
}


// class Glacier
class Glacier
{
 public:
  Glacier();
  ~Glacier();
  void SetGlacierSurface(GlacierSurface *theGlacierSurface) { ptrSurface = theGlacierSurface; }
  GlacierSurface *GetGlacierSurface() const { return ptrSurface; }
  void SetSnow(Snow *theSnow) { ptrSnow = theSnow; }
  Snow *GetSnow() const { return ptrSnow; }
  //  void SetGlacierIce(GlacierIce *theGlacierIce) { ptrIce = theGlacierIce; }
  //  GlacierIce *GetGlacierIce() const { return ptrIce; }
  void SetHBV(HBV *theHBV) { ptrHbv = theHBV; }
  HBV *GetHBV() const { return ptrHbv; }
  void SetKinematicWave(KinematicWave *theKinematicWave) { ptrKwa = theKinematicWave; }
  KinematicWave *GetKinematicWave() const { return ptrKwa; }
  void SetAreaFraction(double value) { areaFraction = value; }
  double GetAreaFraction() const { return areaFraction; }
  void SetSnowStore(double value);
  void SetSubSurfaceHbvStore(double sm, double uz, double lz);
  //  void SetInitialKwaValues();
  void WaterBalance(int timeStep, int initialTimeSteps, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge);
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetSnowCoverFraction() const;
  double GetSnowStore() const;
  double GetMeltWater() const;
  double GetGlacierMassBalance() const;
  double GetKiWaSoilMoisture(double lengthFraction) const;
  double GetKiWaGroundWaterDepth(double lengthFraction) const;
  double GetHbvSoilMoisture() const;
  double GetHbvSoilMoistureDeficit() const;
  double GetHbvPercSoilUpper() const;
  double GetHbvUpperZone() const;
  double GetHbvLowerZone() const;
  double GetLowerRunoff() const;
  double GetUpperRunoff() const;
  double GetRunoff() const;

 private:
  GlacierSurface *ptrSurface;
  Snow *ptrSnow;
  //  GlacierIce *ptrIce;
  HBV *ptrHbv;
  KinematicWave *ptrKwa;
  double areaFraction;
  double glacierMassBalance;               /*  Glacier mass balance (m/timestep)  */
};

Glacier::Glacier():
  glacierMassBalance(0.0)
{
  SetAreaFraction(0.0);
  SetKinematicWave(0);
}

Glacier::~Glacier()
{
}

void Glacier::SetSnowStore(double value) 
{
  GetSnow()->SetSnowStore(value);
}

void Glacier::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
  GetHBV()->SetSubSurfaceHbvStore(sm, uz, lz);
}

/*void Glacier::SetInitialKwaValues()
{
  GetKinematicWave()->SetInitialKwaValues();
}*/

void Glacier::WaterBalance(int timeStep, int initialTimeSteps, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge)
{
  GetGlacierSurface()->WaterBalance(timeStep, initialTimeSteps);
  GetSnow()->WaterBalance(timeStep, GetGlacierSurface()->GetPrecipitation());
  //  GetGlacierIce()->WaterBalance(timeStep, GetGlacierSurface()->GetWaterOutput());
  // Algorithm to be performed in case: no input to landscape element from upstream elements
  /*  GetHBV()->WaterBalance(timeStep, GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction()) 
      + GetSnow()->GetWaterOutput(), missingData);*/
  // Algorithm to be performed in case: input to landscape element from upstream elements
  if (timeStep > initialTimeSteps)
  {
    glacierMassBalance = glacierMassBalance + GetSnow()->GetSnowWaterEquivalentChange() - GetGlacierSurface()->GetIceMelt() * (1.0 - GetSnow()->GetSnowCoverFraction());
  }
  else
  {
    glacierMassBalance = 0.0;
  }
  if (GetHBV())
    GetHBV()->WaterBalance(timeStep, GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction()) 
                           + GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), missingData, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
  if (GetKinematicWave())
    GetKinematicWave()->WaterBalance(timeStep, GetGlacierSurface()->GetIceMelt()*(1.0-GetSnow()->GetSnowCoverFraction()) 
                           + GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), missingData, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
}

double Glacier::GetPrecipitation() const
{
  double precipitation=GetGlacierSurface()->GetPrecipitation();
  return precipitation*GetAreaFraction()/100.0;
}

double Glacier::GetTemperature() const
{
  double temperature=GetGlacierSurface()->GetTemperature();
  return temperature*GetAreaFraction()/100.0;
}

double Glacier::GetSnowCoverFraction() const
{
  double snowCoverFraction=GetSnow()->GetSnowCoverFraction();
  return snowCoverFraction*GetAreaFraction()/100.0;
}

double Glacier::GetSnowStore() const
{
  double snowStore=GetSnow()->GetSnowStore();
  return snowStore*GetAreaFraction()/100.0;
}

double Glacier::GetMeltWater() const
{
  double meltWater=GetSnow()->GetMeltWater();
  return meltWater*GetAreaFraction()/100.0;
}

double Glacier::GetGlacierMassBalance() const
{
    return glacierMassBalance*GetAreaFraction()/100.0;
}

double Glacier::GetKiWaSoilMoisture(double lengthFraction) const
{
  double soilMoisture=GetKinematicWave()->GetSoilMoisture(lengthFraction);
  return soilMoisture*GetAreaFraction()/100.0;
}

double Glacier::GetKiWaGroundWaterDepth(double lengthFraction) const
{
  double groundWaterDepth=GetKinematicWave()->GetGroundWaterDepth(lengthFraction);
  return groundWaterDepth*GetAreaFraction()/100.0;
}

double Glacier::GetHbvSoilMoisture() const
{
  double soilMoisture=GetHBV()->GetSoilMoisture();
  return soilMoisture*GetAreaFraction()/100.0;
}

double Glacier::GetHbvSoilMoistureDeficit() const
{
  double soilMoistureDeficit=GetHBV()->GetSoilMoistureDeficit();
  return soilMoistureDeficit*GetAreaFraction()/100.0;
}

double Glacier::GetHbvPercSoilUpper() const
{
  double percSoilUpper=GetHBV()->GetPercSoilUpper();
  return percSoilUpper*GetAreaFraction()/100.0;
}

double Glacier::GetHbvUpperZone() const
{
  double upperZone=GetHBV()->GetUpperZone();
  return upperZone*GetAreaFraction()/100.0;
}

double Glacier::GetHbvLowerZone() const
{
  double lowerZone=GetHBV()->GetLowerZone();
  return lowerZone*GetAreaFraction()/100.0;
}

double Glacier::GetRunoff() const
{
  double runoff=0.0;
  if (GetHBV()) runoff=runoff+GetHBV()->GetRunoff();
  if (GetKinematicWave()) runoff=runoff+GetKinematicWave()->GetRunoff();
  //  cout << "glacier runoff " << runoff << "  " << GetGlacierSurface()->GetIceMelt() + GetSnow()->GetWaterOutput() << endl;
  return runoff*GetAreaFraction()/100.0;
}

double Glacier::GetLowerRunoff() const
{
  double lowerRunoff=0.0;
  if (GetHBV()) lowerRunoff=lowerRunoff+GetHBV()->GetLowerRunoff();
  if (GetKinematicWave()) lowerRunoff=lowerRunoff+GetKinematicWave()->GetLowerRunoff();
  return lowerRunoff*GetAreaFraction()/100.0;
}

double Glacier::GetUpperRunoff() const
{
  double upperRunoff=0.0;
  if (GetHBV()) upperRunoff=upperRunoff+GetHBV()->GetUpperRunoff();
  if (GetKinematicWave()) upperRunoff=upperRunoff+GetKinematicWave()->GetUpperRunoff();
  return upperRunoff*GetAreaFraction()/100.0;
}


// class HbvAquifer
class HbvAquifer
{
 public:
  HbvAquifer();
  ~HbvAquifer();
  void SetNextHbvAquifer(HbvAquifer *theHbvAquifer) { nextHbvAquifer = theHbvAquifer; }
  HbvAquifer *GetNextHbvAquifer() const { return nextHbvAquifer; }
  void SetVegetation(Vegetation *theVegetation) { ptrVeg = theVegetation; }
  Vegetation *GetVegetation() const { return ptrVeg; }
  void SetSnow(Snow *theSnow) { ptrSnow = theSnow; }
  Snow *GetSnow() const { return ptrSnow; }
  void SetHBV(HBV *theHBV) { ptrHbv = theHBV; }
  HBV *GetHBV() const { return ptrHbv; }
  void SetAreaFraction(double value) { areaFraction = value; }
  double GetAreaFraction() const { return areaFraction; }
  void SetSnowStore(double value);
  void SetSubSurfaceHbvStore(double sm, double uz, double lz);
  // Algorithm to be performed in case: no input to landscape element from upstream elements
  //  void WaterBalance(int timeStep);
  // Algorithm to be performed in case: input to landscape element from upstream elements
  void WaterBalance(int timeStep, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge);
  double GetTotalHbvAreaFraction() const;
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetSnowCoverFraction() const;
  double GetSnowStore() const;
  double GetMeltWater() const;
  double GetSoilMoisture() const;
  double GetSoilMoistureDeficit() const;
  double GetPercSoilUpper() const;
  double GetUpperZone() const;
  double GetLowerZone() const;
  double GetInterceptionLoss() const;
  double GetTranspSoilEvap() const;
  double GetLowerRunoff() const;
  double GetUpperRunoff() const;
  double GetRunoff() const;

 private:
  Vegetation *ptrVeg;
  Snow *ptrSnow;
  HBV *ptrHbv;
  HbvAquifer *nextHbvAquifer;
  double areaFraction;
};

HbvAquifer::HbvAquifer()
{
  SetAreaFraction(0.0);
  SetNextHbvAquifer(0);
}

HbvAquifer::~HbvAquifer()
{
}

void HbvAquifer::SetSnowStore(double value)
{
  if (GetNextHbvAquifer()) 
    GetNextHbvAquifer()->SetSnowStore(value);
  GetSnow()->SetSnowStore(value);
}

void HbvAquifer::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
  if (GetNextHbvAquifer()) 
    GetNextHbvAquifer()->SetSubSurfaceHbvStore(sm, uz, lz);
  GetHBV()->SetSubSurfaceHbvStore(sm, uz, lz);
}

// Algorithm to be performed in case: no input to landscape element from upstream elements
/*void HbvAquifer::WaterBalance(int timeStep)
  {
  GetVegetation()->WaterBalance(timeStep);
  GetSnow()->WaterBalance(timeStep, GetVegetation()->GetThroughFall());
  GetHBV()->WaterBalance(timeStep, GetSnow()->GetWaterOutput(), GetVegetation()->GetDryPeriod());
  }*/
// End algorithm to be performed in case: no input to landscape element from upstream elements

// Algorithm to be performed in case: input to landscape element from upstream elements
void HbvAquifer::WaterBalance(int timeStep, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge)
{
  GetVegetation()->WaterBalance(timeStep);
  GetSnow()->WaterBalance(timeStep, GetVegetation()->GetThroughFall());
  GetHBV()->WaterBalance(timeStep, GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), GetVegetation()->GetDryPeriod(), upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
}
// End algorithm to be performed in case: input to landscape element from upstream elements

double HbvAquifer::GetTotalHbvAreaFraction() const
{
  double totalHbvAreaFraction=0.0;
  if (GetNextHbvAquifer()) 
    totalHbvAreaFraction=totalHbvAreaFraction+GetNextHbvAquifer()->GetTotalHbvAreaFraction();
  totalHbvAreaFraction=totalHbvAreaFraction+GetAreaFraction();
  return totalHbvAreaFraction;
}

double HbvAquifer::GetPrecipitation() const
{
  double precipitation=0.0;
  if (GetNextHbvAquifer()) 
    precipitation=precipitation+GetNextHbvAquifer()->GetPrecipitation();
  precipitation=precipitation+(GetVegetation()->GetPrecipitation()*GetAreaFraction()/100.0);
  return precipitation;
}

double HbvAquifer::GetTemperature() const
{
  double temperature=0.0;
  if (GetNextHbvAquifer()) 
    temperature=temperature+GetNextHbvAquifer()->GetTemperature();
  temperature=temperature+(GetVegetation()->GetTemperature()*GetAreaFraction()/100.0);
  return temperature;
}

double HbvAquifer::GetSnowCoverFraction() const
{
  double snowCoverFraction=0.0;
  if (GetNextHbvAquifer()) 
    snowCoverFraction=snowCoverFraction+GetNextHbvAquifer()->GetSnowCoverFraction();
  snowCoverFraction=snowCoverFraction+(GetSnow()->GetSnowCoverFraction()*GetAreaFraction()/100.0);
  return snowCoverFraction;
}

double HbvAquifer::GetSnowStore() const
{
  double snowStore=0.0;
  if (GetNextHbvAquifer()) 
    snowStore=snowStore+GetNextHbvAquifer()->GetSnowStore();
  snowStore=snowStore+(GetSnow()->GetSnowStore()*GetAreaFraction()/100.0);
  return snowStore;
}

double HbvAquifer::GetMeltWater() const
{
  double meltWater=0.0;
  if (GetNextHbvAquifer()) 
    meltWater=meltWater+GetNextHbvAquifer()->GetMeltWater();
  meltWater=meltWater+(GetSnow()->GetMeltWater()*GetAreaFraction()/100.0);
  return meltWater;
}

double HbvAquifer::GetSoilMoisture() const
{
  double soilMoisture=0.0;
  if (GetNextHbvAquifer()) 
    soilMoisture=soilMoisture+GetNextHbvAquifer()->GetSoilMoisture();
  soilMoisture=soilMoisture+(GetHBV()->GetSoilMoisture()*GetAreaFraction()/100.0);
  return soilMoisture;
}

double HbvAquifer::GetSoilMoistureDeficit() const
{
  double soilMoistureDeficit=0.0;
  if (GetNextHbvAquifer()) 
    soilMoistureDeficit=soilMoistureDeficit+GetNextHbvAquifer()->GetSoilMoistureDeficit();
  soilMoistureDeficit=soilMoistureDeficit+(GetHBV()->GetSoilMoistureDeficit()*GetAreaFraction()/100.0);
  return soilMoistureDeficit;
}

double HbvAquifer::GetPercSoilUpper() const
{
  double percSoilUpper=0.0;
  if (GetNextHbvAquifer()) 
    percSoilUpper=percSoilUpper+GetNextHbvAquifer()->GetPercSoilUpper();
  percSoilUpper=percSoilUpper+(GetHBV()->GetPercSoilUpper()*GetAreaFraction()/100.0);
  return percSoilUpper;
}

double HbvAquifer::GetUpperZone() const
{
  double upperZone=0.0;
  if (GetNextHbvAquifer()) 
    upperZone=upperZone+GetNextHbvAquifer()->GetUpperZone();
  upperZone=upperZone+(GetHBV()->GetUpperZone()*GetAreaFraction()/100.0);
  return upperZone;
}

double HbvAquifer::GetLowerZone() const
{
  double lowerZone=0.0;
  if (GetNextHbvAquifer()) 
    lowerZone=lowerZone+GetNextHbvAquifer()->GetLowerZone();
  lowerZone=lowerZone+(GetHBV()->GetLowerZone()*GetAreaFraction()/100.0);
  return lowerZone;
}

double HbvAquifer::GetInterceptionLoss() const
{
  double interceptionLoss=0.0;
  if (GetNextHbvAquifer()) 
    interceptionLoss=interceptionLoss+GetNextHbvAquifer()->GetInterceptionLoss();
  interceptionLoss=interceptionLoss+(GetVegetation()->GetInterceptionLoss()*GetAreaFraction()/100.0);
  return interceptionLoss;
}

double HbvAquifer::GetTranspSoilEvap() const
{
  double transpSoilEvap=0.0;
  if (GetNextHbvAquifer()) 
    transpSoilEvap=transpSoilEvap+GetNextHbvAquifer()->GetTranspSoilEvap();
  transpSoilEvap=transpSoilEvap+(GetHBV()->GetTranspSoilEvap()*GetAreaFraction()/100.0);
  return transpSoilEvap;
}

double HbvAquifer::GetRunoff() const
{
  double runoff=0.0;
  if (GetNextHbvAquifer()) 
    runoff=runoff+GetNextHbvAquifer()->GetRunoff();
  runoff=runoff+(GetHBV()->GetRunoff()*GetAreaFraction()/100.0);
  return runoff;
}

double HbvAquifer::GetLowerRunoff() const
{
  double lowerRunoff=0.0;
  if (GetNextHbvAquifer()) 
    lowerRunoff=lowerRunoff+GetNextHbvAquifer()->GetLowerRunoff();
  lowerRunoff=lowerRunoff+(GetHBV()->GetLowerRunoff()*GetAreaFraction()/100.0);
  return lowerRunoff;
}

double HbvAquifer::GetUpperRunoff() const
{
  double upperRunoff=0.0;
  if (GetNextHbvAquifer()) 
    upperRunoff=upperRunoff+GetNextHbvAquifer()->GetUpperRunoff();
  upperRunoff=upperRunoff+(GetHBV()->GetUpperRunoff()*GetAreaFraction()/100.0);
  return upperRunoff;
}


// class KwaAquifer
class KwaAquifer
{
 public:
  KwaAquifer();
  ~KwaAquifer();
  void SetNextKwaAquifer(KwaAquifer *theKwaAquifer) { nextKwaAquifer = theKwaAquifer; }
  KwaAquifer *GetNextKwaAquifer() const { return nextKwaAquifer; }
  void SetVegetation(Vegetation *theVegetation) { ptrVeg = theVegetation; }
  Vegetation *GetVegetation() const { return ptrVeg; }
  void SetSnow(Snow *theSnow) { ptrSnow = theSnow; }
  Snow *GetSnow() const { return ptrSnow; }
  void SetKinematicWave(KinematicWave *theKinematicWave) { ptrKwa = theKinematicWave; }
  KinematicWave *GetKinematicWave() const { return ptrKwa; }
  void SetAreaFraction(double value) { areaFraction = value; }
  double GetAreaFraction() const { return areaFraction; }
  void SetSnowStore(double value);
  void WaterBalance(int timeStep, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge);
  double GetTotalKwaAreaFraction() const;
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetSnowCoverFraction() const;
  double GetSnowStore() const;
  double GetMeltWater() const;
  double GetSoilMoisture(double lengthFraction) const;
  double GetGroundWaterDepth(double lengthFraction) const;
  double GetInterceptionLoss() const;
  double GetTranspSoilEvap() const;
  double GetLowerRunoff() const;
  double GetUpperRunoff() const;
  double GetRunoff() const;

 private:
  Vegetation *ptrVeg;
  Snow *ptrSnow;
  KinematicWave *ptrKwa;
  KwaAquifer *nextKwaAquifer;
  double areaFraction;
};

KwaAquifer::KwaAquifer()
{
  SetAreaFraction(0);
  SetNextKwaAquifer(0);
  SetKinematicWave(0);
}

KwaAquifer::~KwaAquifer()
{
}

void KwaAquifer::SetSnowStore(double value)
{
  if (GetNextKwaAquifer()) 
    GetNextKwaAquifer()->SetSnowStore(value);
  GetSnow()->SetSnowStore(value);
}

void KwaAquifer::WaterBalance(int timeStep, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge)
{
  //  cout << "kwaAquifer\n";
  GetVegetation()->WaterBalance(timeStep);
  GetSnow()->WaterBalance(timeStep, GetVegetation()->GetThroughFall());
  GetKinematicWave()->WaterBalance(timeStep, GetSnow()->GetWaterOutput(), GetSnow()->GetSnowCoverFraction(), GetVegetation()->GetDryPeriod(), upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
}

double KwaAquifer::GetTotalKwaAreaFraction() const
{
  double totalKwaAreaFraction=0.0;
  if (GetNextKwaAquifer()) 
    totalKwaAreaFraction=totalKwaAreaFraction+GetNextKwaAquifer()->GetTotalKwaAreaFraction();
  totalKwaAreaFraction=totalKwaAreaFraction+GetAreaFraction();
  return totalKwaAreaFraction;
}

double KwaAquifer::GetPrecipitation() const
{
  double precipitation=0.0;
  if (GetNextKwaAquifer()) 
    precipitation=precipitation+GetNextKwaAquifer()->GetPrecipitation();
  precipitation=precipitation+(GetVegetation()->GetPrecipitation()*GetAreaFraction()/100.0);
  return precipitation;
}

double KwaAquifer::GetTemperature() const
{
  double temperature=0.0;
  if (GetNextKwaAquifer()) 
    temperature=temperature+GetNextKwaAquifer()->GetTemperature();
  temperature=temperature+(GetVegetation()->GetTemperature()*GetAreaFraction()/100.0);
  return temperature;
}

double KwaAquifer::GetSnowCoverFraction() const
{
  double snowCoverFraction=0.0;
  if (GetNextKwaAquifer()) 
    snowCoverFraction=snowCoverFraction+GetNextKwaAquifer()->GetSnowCoverFraction();
  snowCoverFraction=snowCoverFraction+(GetSnow()->GetSnowCoverFraction()*GetAreaFraction()/100.0);
  return snowCoverFraction;
}

double KwaAquifer::GetSnowStore() const
{
  double snowStore=0.0;
  if (GetNextKwaAquifer()) 
    snowStore=snowStore+GetNextKwaAquifer()->GetSnowStore();
  snowStore=snowStore+(GetSnow()->GetSnowStore()*GetAreaFraction()/100.0);
  return snowStore;
}

double KwaAquifer::GetMeltWater() const
{
  double meltWater=0.0;
  if (GetNextKwaAquifer()) 
    meltWater=meltWater+GetNextKwaAquifer()->GetMeltWater();
  meltWater=meltWater+(GetSnow()->GetMeltWater()*GetAreaFraction()/100.0);
  return meltWater;
}

double KwaAquifer::GetSoilMoisture(double lengthFraction) const
{
  double soilMoisture=0.0;
  if (GetNextKwaAquifer()) 
    soilMoisture=soilMoisture+GetNextKwaAquifer()->GetSoilMoisture(lengthFraction);
  soilMoisture=soilMoisture+(GetKinematicWave()->GetSoilMoisture(lengthFraction)*GetAreaFraction()/100.0);
  return soilMoisture;
}

double KwaAquifer::GetGroundWaterDepth(double lengthFraction) const
{
  double groundWaterDepth=0.0;
  if (GetNextKwaAquifer()) 
    groundWaterDepth=groundWaterDepth+GetNextKwaAquifer()->GetGroundWaterDepth(lengthFraction);
  groundWaterDepth=groundWaterDepth+(GetKinematicWave()->GetGroundWaterDepth(lengthFraction)*GetAreaFraction()/100.0);
  return groundWaterDepth;
}

double KwaAquifer::GetInterceptionLoss() const
{
  double interceptionLoss=0.0;
  if (GetNextKwaAquifer()) 
    interceptionLoss=interceptionLoss+GetNextKwaAquifer()->GetInterceptionLoss();
  interceptionLoss=interceptionLoss+(GetVegetation()->GetInterceptionLoss()*GetAreaFraction()/100.0);
  return interceptionLoss;
}

double KwaAquifer::GetTranspSoilEvap() const
{
  double transpSoilEvap=0.0;
  if (GetNextKwaAquifer()) 
    transpSoilEvap=transpSoilEvap+GetNextKwaAquifer()->GetTranspSoilEvap();
  transpSoilEvap=transpSoilEvap+(GetKinematicWave()->GetTranspSoilEvap()*GetAreaFraction()/100.0);
  return transpSoilEvap;
}

double KwaAquifer::GetRunoff() const
{
  double runoff=0;
  if (GetNextKwaAquifer()) 
    runoff=runoff+GetNextKwaAquifer()->GetRunoff();
  runoff=runoff+GetKinematicWave()->GetRunoff()*GetAreaFraction()/100.0;
  return runoff;
}

double KwaAquifer::GetLowerRunoff() const
{
  double lowerRunoff=0;
  if (GetNextKwaAquifer()) 
    lowerRunoff=lowerRunoff+GetNextKwaAquifer()->GetLowerRunoff();
  lowerRunoff=lowerRunoff+GetKinematicWave()->GetLowerRunoff()*GetAreaFraction()/100.0;
  return lowerRunoff;
}

double KwaAquifer::GetUpperRunoff() const
{
  double upperRunoff=0;
  if (GetNextKwaAquifer()) 
    upperRunoff=upperRunoff+GetNextKwaAquifer()->GetUpperRunoff();
  upperRunoff=upperRunoff+GetKinematicWave()->GetUpperRunoff()*GetAreaFraction()/100.0;
  return upperRunoff;
}


// class DistributedElement
class DistributedElement
{
 public:
  DistributedElement();
  ~DistributedElement(); 
  void SetGeoIndex(int value) { geoIndex = value; }
  int GetGeoIndex() const { return geoIndex; }
  void SetLandIndex(int value) { landIndex = value; }
  int GetLandIndex() const { return landIndex; }
  void SetLakeNumber(int value) { lakeNumber = value; }
  int GetLakeNumber() const { return lakeNumber; }
  void AllocateInputValues(int value);
  int GetNumberInputValues() const { return numberInputValues; }
  void SetInputValue(int i, double value) { inputArray[i] = value; }
  double GetInputValue(int i) const { return inputArray[i]; }
  void SetArea(double value) { area = value; }
  double GetArea() const { return area; }
  void SetElevation(double value) { elevation = value; }
  double GetElevation() const { return elevation; }
  void SetSlopeLength(double value) { slopeLength = value; }
  double GetSlopeLength() const { return slopeLength; }
  void SetSlopeAngle(double value);
  double GetSlopeAngle() const { return slopeAngle; }
  void SetAspect(double value) { aspect = value; }
  double GetAspect() const { return aspect; }
  void SetFlowDirection(int value) { flowDirection = value; }
  int GetFlowDirection() const { return flowDirection; }
  void SetSelectedKiWaTimeSeriesElements(SelectedKiWaTimeSeriesElements *object) { selectedKiWaElements = object; }
  SelectedKiWaTimeSeriesElements *GetSelectedKiWaTimeSeriesElements() const { return selectedKiWaElements; }
  void SetSelectedHbvTimeSeriesElements(SelectedHbvTimeSeriesElements *object) { selectedHbvElements = object; }
  SelectedHbvTimeSeriesElements *GetSelectedHbvTimeSeriesElements() const { return selectedHbvElements; }
  void SetCommonPar(ParametersCommon *parObj) { commonPar = parObj; }
  ParametersCommon *GetCommonPar() const { return commonPar; }
  void SetPrecipitationCorrection(double value) { precipitationCorrection = value; }
  void SetTotalReservoir(TotalReservoirStorage *resObj) { totalReservoir = resObj; }
  TotalReservoirStorage *GetTotalReservoir() const { return totalReservoir; }
  double GetPrecipitationCorrection() const { return precipitationCorrection; }
  void SetTemperatureCorrection(double value) { temperatureCorrection = value; }
  double GetTemperatureCorrection() const { return temperatureCorrection; }
  void SetNumUpLand(int value);
  int GetNumUpLand() const { return numUpLand; }
  void SetUpLandFlow(int k, DistributedElement *theElement) {upLandFlow[k] = theElement; }
  DistributedElement *GetUpLandFlow(int k) const { return upLandFlow[k]; }
  void SetNextElement(DistributedElement *theElement) { nextElement = theElement; }
  DistributedElement *GetNextElement() const { return nextElement; }
  void SetWaterCourseElement(WaterCourse *theElement) { waterCourseElement = theElement; }
  WaterCourse *GetWaterCourseElement() const { return waterCourseElement; }
  void SetLake(Lake *theLake) { lake = theLake; }
  Lake *GetLake() const { return lake; }
  void SetGlacier(Glacier *theGlacier) { glacier = theGlacier; }
  Glacier *GetGlacier() const { return glacier; }
  void SetHbvAquifer(HbvAquifer *theHbvAquifer) { hbvAquifer = theHbvAquifer; }
  HbvAquifer *GetHbvAquifer() const { return hbvAquifer; }
  void SetKwaAquifer(KwaAquifer *theKwaAquifer) { kwaAquifer = theKwaAquifer; }
  KwaAquifer *GetKwaAquifer() const { return kwaAquifer; }
  void SetPrecStationsWeightedElevation(double value) { precStationsWeightedElevation = value; }
  double GetPrecStationsWeightedElevation() { return precStationsWeightedElevation; }
  void SetTempStationsWeightedElevation(double value) { tempStationsWeightedElevation = value; }
  double GetTempStationsWeightedElevation() { return tempStationsWeightedElevation; }
  void AllocateMetSeries(int numberPrecSeries, int numberTempSeries);
  void SetMetSeriesNumber(int k, int number) { metSeriesNumber[k] = number; }
  void SetMetSeriesWeight(int k, double weight) { metSeriesWeight[k] = weight; }
  int GetMetSeriesNumber(int k) { return metSeriesNumber[k]; }
  double GetMetSeriesWeight(int k) { return metSeriesWeight[k]; }
  void AllocateKiWaWaterBalance(int numberTimeSteps); 
  void AllocateHbvWaterBalance(int numberTimeSteps); 
  void SetDistributedElementPrecipitation(int index, double value) { distributedElementPrecipitation[index] = value; }
  double GetDistributedElementPrecipitation(int index) const { return distributedElementPrecipitation[index]; }
  void SetDistributedElementTemperature(int index, double value) { distributedElementTemperature[index] = value; }
  double GetDistributedElementTemperature(int index) const { return distributedElementTemperature[index]; }
  void SetDistributedElementSnowStore(int index, double value) { distributedElementSnowStore[index] = value; }
  double GetDistributedElementSnowStore(int index) const { return distributedElementSnowStore[index]; }
  void SetDistributedElementSnowCoverFraction(int index, double value) { distributedElementSnowCoverFraction[index] = value; }
  double GetDistributedElementSnowCoverFraction(int index) const { return distributedElementSnowCoverFraction[index]; }
  void SetDistributedElementMeltWater(int index, double value) { distributedElementMeltWater[index] = value; }
  double GetDistributedElementMeltWater(int index) const { return distributedElementMeltWater[index]; }
  void SetDistributedElementGlacierMassBalance(int index, double value) { distributedElementGlacierMassBalance[index] = value; }
  double GetDistributedElementGlacierMassBalance(int index) const { return distributedElementGlacierMassBalance[index]; }
  void SetDistributedElementEvapotranspiration(int index, double value) { distributedElementEvapotranspiration[index] = value; }
  double GetDistributedElementEvapotranspiration(int index) const { return distributedElementEvapotranspiration[index]; }
  void SetDistributedElementRunoff(int index, double value) { distributedElementRunoff[index] = value; }
  double GetDistributedElementRunoff(int index) const { return distributedElementRunoff[index]; }
  void SetDistributedElementKiWaSoilMoistureOne(int index, double value) { distributedElementKiWaSoilMoistureOne[index] = value; }
  double GetDistributedElementKiWaSoilMoistureOne(int index) const { return distributedElementKiWaSoilMoistureOne[index]; }
  void SetDistributedElementKiWaSoilMoistureTwo(int index, double value) { distributedElementKiWaSoilMoistureTwo[index] = value; }
  double GetDistributedElementKiWaSoilMoistureTwo(int index) const { return distributedElementKiWaSoilMoistureTwo[index]; }
  void SetDistributedElementKiWaGroundWaterDepthOne(int index, double value) { distributedElementKiWaGroundWaterDepthOne[index] = value; }
  double GetDistributedElementKiWaGroundWaterDepthOne(int index) const { return distributedElementKiWaGroundWaterDepthOne[index]; }
  void SetDistributedElementKiWaGroundWaterDepthTwo(int index, double value) { distributedElementKiWaGroundWaterDepthTwo[index] = value; }
  double GetDistributedElementKiWaGroundWaterDepthTwo(int index) const { return distributedElementKiWaGroundWaterDepthTwo[index]; }
  void SetDistributedElementHbvSoilMoisture(int index, double value) { distributedElementHbvSoilMoisture[index] = value; }
  double GetDistributedElementHbvSoilMoisture(int index) const { return distributedElementHbvSoilMoisture[index]; }
  void SetDistributedElementHbvSoilMoistureDeficit(int index, double value)  
    { distributedElementHbvSoilMoistureDeficit[index] = value; }
  double GetDistributedElementHbvSoilMoistureDeficit(int index) const { return distributedElementHbvSoilMoistureDeficit[index]; }
  void SetDistributedElementHbvPercSoilUpper(int index, double value) { distributedElementHbvPercSoilUpper[index] = value; }
  double GetDistributedElementHbvPercSoilUpper(int index) const { return distributedElementHbvPercSoilUpper[index]; }
  void SetDistributedElementHbvUpperZone(int index, double value) { distributedElementHbvUpperZone[index] = value; }
  double GetDistributedElementHbvUpperZone(int index) const { return distributedElementHbvUpperZone[index]; }
  void SetDistributedElementHbvLowerZone(int index, double value) { distributedElementHbvLowerZone[index] = value; }
  double GetDistributedElementHbvLowerZone(int index) const { return distributedElementHbvLowerZone[index]; }
  void SetAccumulatedDischarge(double upLandValue, double localValue, bool flowHierarchy, bool forceDirect);
  double GetAccumulatedDischarge() const { return accumulatedDischarge; }
  void SetAccumulatedLowerDischarge(double upLandValue, double localValue, bool flowHierarchy, bool forceDirect);
  double GetAccumulatedLowerDischarge() const { return accumulatedLowerDischarge; }
  void SetAccumulatedUpperDischarge(double upLandValue, double localValue, bool flowHierarchy, bool forceDirect);
  double GetAccumulatedUpperDischarge() const { return accumulatedUpperDischarge; }
  void SetAccumulatedSum(double value) { accumulatedSum = value; }
  double GetAccumulatedSum() const { return accumulatedSum; }
  void SetAccumulatedSumSnow(double value) { accumulatedSumSnow = value; }
  double GetAccumulatedSumSnow() const { return accumulatedSumSnow; }
  void SetAccumulatedSumGlacier(double value) { accumulatedSumGlacier = value; }
  double GetAccumulatedSumGlacier() const { return accumulatedSumGlacier; }
  void SetAccumulatedSumHbv(double value) { accumulatedSumHbv = value; }
  double GetAccumulatedSumHbv() const { return accumulatedSumHbv; }
  void SetAccumulatedPrecipitation(double value) { accumulatedPrecipitation = value; }
  double GetAccumulatedPrecipitation() const { return accumulatedPrecipitation; }
  void SetAccumulatedTemperature(double value) { accumulatedTemperature = value; }
  double GetAccumulatedTemperature() const { return accumulatedTemperature; }
  void SetAccumulatedSnowStore(double value) { accumulatedSnowStore = value; }
  double GetAccumulatedSnowStore() const { return accumulatedSnowStore; }
  void SetAccumulatedGlacierMassBalance(double value) { accumulatedGlacierMassBalance = value; }
  double GetAccumulatedGlacierMassBalance() const { return accumulatedGlacierMassBalance; }
  void SetAccumulatedEvapotranspiration(double value) { accumulatedEvapotranspiration = value; }
  double GetAccumulatedEvapotranspiration() const { return accumulatedEvapotranspiration; }
  void SetAccumulatedRunoff(double value) { accumulatedRunoff = value; }
  double GetAccumulatedRunoff() const { return accumulatedRunoff; }
  void SetAccumulatedHbvSoilMoisture(double value) { accumulatedHbvSoilMoisture = value; }
  double GetAccumulatedHbvSoilMoisture() const { return accumulatedHbvSoilMoisture; }
  void SetAccumulatedHbvSoilMoistureDeficit(double value) { accumulatedHbvSoilMoistureDeficit = value; }
  double GetAccumulatedHbvSoilMoistureDeficit() const { return accumulatedHbvSoilMoistureDeficit; }
  void SetAccumulatedHbvPercSoilUpper(double value) { accumulatedHbvPercSoilUpper = value; }
  double GetAccumulatedHbvPercSoilUpper() const { return accumulatedHbvPercSoilUpper; }
  void SetAccumulatedHbvUpperZone(double value) { accumulatedHbvUpperZone = value; }
  double GetAccumulatedHbvUpperZone() const { return accumulatedHbvUpperZone; }
  void SetAccumulatedHbvLowerZone(double value) { accumulatedHbvLowerZone = value; }
  double GetAccumulatedHbvLowerZone() const { return accumulatedHbvLowerZone; }
  void SetSnowStore(double value);
  void SetSubSurfaceHbvStore(double sm, double uz, double lz);
  double GetLandArea() const;
  // Algorithm to be performed in case: no input to landscape element from upstream elements
  //  void WaterBalance(int timeStep) const;
  // Algorithm to be performed in case: input to landscape element from upstream elements
  void WaterBalance(int timeStep, int initialTimeSteps, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge) const;
  double GetPrecipitation() const;
  double GetTemperature() const;
  double GetSnowCoverFraction() const;
  double GetSnowStore() const;
  double GetMeltWater() const;
  double GetGlacierMassBalance() const;
  double GetKiWaSoilMoisture(double lengthFraction) const;
  double GetKiWaGroundWaterDepth(double lengthFraction) const;
  double GetHbvSoilMoisture() const;
  double GetHbvSoilMoistureDeficit() const;
  double GetHbvPercSoilUpper() const;
  double GetHbvUpperZone() const;
  double GetHbvLowerZone() const;
  double GetInterceptionLoss() const;
  double GetTranspSoilEvap() const;
  double GetLakeEvap() const;
  double GetRunoff() const;
  double GetLowerRunoff() const;
  double GetUpperRunoff() const;
  double GetDischarge() const;
  double GetLowerDischarge() const;
  double GetUpperDischarge() const;
  void SetTotalReservoirStorage(DateTimeInformation * const DateTimeStore, int index, bool * firstTotal, bool inputDataFound);

private:
  int geoIndex;
  int landIndex;
  int lakeNumber;
  int numUpLand;
  int numberInputValues;
  int flowDirection;
  double * inputArray;
  double accumulatedSum;
  double accumulatedSumSnow;
  double accumulatedSumGlacier;
  double accumulatedSumHbv;
  double area;
  double elevation;
  double slopeLength;
  double slopeAngle;
  double aspect;
  double precipitationCorrection;
  double temperatureCorrection;
  double accumulatedDischarge;
  double accumulatedLowerDischarge;
  double accumulatedUpperDischarge;
  double accumulatedPrecipitation;
  double accumulatedTemperature;
  double accumulatedSnowStore;
  double accumulatedGlacierMassBalance;
  double accumulatedEvapotranspiration;
  double accumulatedRunoff;
  double accumulatedHbvSoilMoisture;
  double accumulatedHbvSoilMoistureDeficit;
  double accumulatedHbvPercSoilUpper;
  double accumulatedHbvUpperZone;
  double accumulatedHbvLowerZone;
  double precStationsWeightedElevation;
  double tempStationsWeightedElevation;
  int * metSeriesNumber;
  double * metSeriesWeight;
  double *distributedElementPrecipitation;
  double *distributedElementTemperature;
  double *distributedElementSnowStore;           
  double *distributedElementSnowCoverFraction;           
  double *distributedElementMeltWater;           
  double *distributedElementGlacierMassBalance;
  double *distributedElementEvapotranspiration;  
  double *distributedElementRunoff;              
  double *distributedElementKiWaSoilMoistureOne;     
  double *distributedElementKiWaSoilMoistureTwo;     
  double *distributedElementKiWaGroundWaterDepthOne;     
  double *distributedElementKiWaGroundWaterDepthTwo;     
  double *distributedElementHbvSoilMoisture;     
  double *distributedElementHbvSoilMoistureDeficit;     
  double *distributedElementHbvPercSoilUpper;     
  double *distributedElementHbvUpperZone;        
  double *distributedElementHbvLowerZone;        
  DistributedElement **upLandFlow;
  DistributedElement *nextElement;
  WaterCourse *waterCourseElement;
  Lake *lake;
  Glacier *glacier;
  HbvAquifer *hbvAquifer;
  KwaAquifer *kwaAquifer;
  SelectedKiWaTimeSeriesElements *selectedKiWaElements;
  SelectedHbvTimeSeriesElements *selectedHbvElements;
  ParametersCommon *commonPar;
  TotalReservoirStorage * totalReservoir;
};

DistributedElement::DistributedElement():
  lakeNumber(-9999),
  numUpLand(0),
  flowDirection(0),
  accumulatedSum(0.0),
  accumulatedSumSnow(0.0),
  accumulatedSumGlacier(0.0),
  accumulatedSumHbv(0.0),
  area(0.0),
  elevation(0.0),
  slopeLength(0.0),
  slopeAngle(0.0),
  aspect(0.0),
  precipitationCorrection(1.0),
  temperatureCorrection(0.0),
  accumulatedDischarge(0.0),
  accumulatedLowerDischarge(0.0),
  accumulatedUpperDischarge(0.0),
  accumulatedPrecipitation(0.0),
  accumulatedTemperature(0.0),
  accumulatedSnowStore(0.0),
  accumulatedGlacierMassBalance(0.0),
  accumulatedEvapotranspiration(0.0),
  accumulatedRunoff(0.0),
  accumulatedHbvSoilMoisture(0.0),
  accumulatedHbvSoilMoistureDeficit(0.0),
  accumulatedHbvPercSoilUpper(0.0),
  accumulatedHbvUpperZone(0.0),
  accumulatedHbvLowerZone(0.0),
  precStationsWeightedElevation(0.0),
  tempStationsWeightedElevation(0.0)
{
  upLandFlow=0;
  SetWaterCourseElement(0); 
  SetNextElement(0);
  SetLake(0);
  SetGlacier(0);
  SetHbvAquifer(0);
  SetKwaAquifer(0);
  SetSelectedKiWaTimeSeriesElements(0);
  SetSelectedHbvTimeSeriesElements(0);
  SetCommonPar(0);
  SetTotalReservoir(0);
}

DistributedElement::~DistributedElement()
{ 
}

void DistributedElement::SetSlopeAngle(double value) 
{
  if (value < 0.01) value=0.01;
  slopeAngle = value;
}

void DistributedElement::SetNumUpLand(int value) 
{
  numUpLand = value;
  DistributedElement **upLand = new DistributedElement * [value];
  upLandFlow = upLand;
}

void DistributedElement::AllocateInputValues(int value)
{
  int i;
  numberInputValues = value;
  inputArray = new double [numberInputValues];
  for (i=0; i<numberInputValues; i++) inputArray[i] = missingData;
}

void DistributedElement::AllocateMetSeries(int numberPrecSeries, int numberTempSeries) 
{
  metSeriesNumber = new int [numberPrecSeries+numberTempSeries];
  metSeriesWeight = new double [numberPrecSeries+numberTempSeries];
}

void DistributedElement::AllocateKiWaWaterBalance(int numberTimeSteps) 
{
  distributedElementPrecipitation = new double [numberTimeSteps];
  distributedElementTemperature = new double [numberTimeSteps];
  distributedElementSnowStore = new double [numberTimeSteps];
  distributedElementSnowCoverFraction = new double [numberTimeSteps];
  distributedElementMeltWater = new double [numberTimeSteps];
  distributedElementGlacierMassBalance = new double[numberTimeSteps];
  distributedElementEvapotranspiration = new double [numberTimeSteps];
  distributedElementRunoff = new double [numberTimeSteps];
  distributedElementKiWaSoilMoistureOne = new double [numberTimeSteps];
  distributedElementKiWaSoilMoistureTwo = new double [numberTimeSteps];
  distributedElementKiWaGroundWaterDepthOne = new double [numberTimeSteps];
  distributedElementKiWaGroundWaterDepthTwo = new double [numberTimeSteps];
}

void DistributedElement::AllocateHbvWaterBalance(int numberTimeSteps) 
{
  distributedElementPrecipitation = new double [numberTimeSteps];
  distributedElementTemperature = new double [numberTimeSteps];
  distributedElementSnowStore = new double [numberTimeSteps];
  distributedElementSnowCoverFraction = new double [numberTimeSteps];
  distributedElementMeltWater = new double [numberTimeSteps];
  distributedElementGlacierMassBalance = new double[numberTimeSteps];
  distributedElementEvapotranspiration = new double [numberTimeSteps];
  distributedElementRunoff = new double [numberTimeSteps];
  distributedElementHbvSoilMoisture = new double [numberTimeSteps];
  distributedElementHbvSoilMoistureDeficit = new double [numberTimeSteps];
  distributedElementHbvPercSoilUpper = new double [numberTimeSteps];
  distributedElementHbvUpperZone = new double [numberTimeSteps];
  distributedElementHbvLowerZone = new double [numberTimeSteps];
}

void DistributedElement::SetAccumulatedDischarge(double upLandValue, double localValue, bool flowHierarchy, bool forceDirect) 
{
  // Algorithm to be performed in case: no input to landscape element from upstream elements
  //  accumulatedDischarge = upLandValue+localValue;
  // Algorithm to be performed in case: input to landscape element from upstream elements
  //  if (flowHierarchy && GetHbvAquifer()) 
  if (flowHierarchy && !forceDirect) 
    accumulatedDischarge = localValue;
  else
    accumulatedDischarge = upLandValue+localValue;
}

void DistributedElement::SetAccumulatedLowerDischarge(double upLandValue, double localValue, bool flowHierarchy, bool forceDirect) 
{
 if (flowHierarchy && !forceDirect) 
    accumulatedLowerDischarge = localValue;
  else
    accumulatedLowerDischarge = upLandValue+localValue;
}

void DistributedElement::SetAccumulatedUpperDischarge(double upLandValue, double localValue, bool flowHierarchy, bool forceDirect) 
{
 if (flowHierarchy && !forceDirect) 
    accumulatedUpperDischarge = localValue;
  else
    accumulatedUpperDischarge = upLandValue+localValue;
}


void DistributedElement::SetSnowStore(double value)
{
  if (GetGlacier()) GetGlacier()->SetSnowStore(value);
  if (GetHbvAquifer()) GetHbvAquifer()->SetSnowStore(value);
  if (GetKwaAquifer()) GetKwaAquifer()->SetSnowStore(value);
}

void DistributedElement::SetSubSurfaceHbvStore(double sm, double uz, double lz)
{
  if (GetGlacier()) 
    if (GetGlacier()->GetHBV()) GetGlacier()->SetSubSurfaceHbvStore(sm, uz, lz);
  if (GetHbvAquifer()) GetHbvAquifer()->SetSubSurfaceHbvStore(sm, uz, lz);
}

double DistributedElement::GetLandArea() const
{
  double sumAreaFraction = 0.0;
  bool elementFound = false;
  if (GetGlacier())
    {
      sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
      elementFound = true;
    }
  if (GetHbvAquifer())
    {
      sumAreaFraction = sumAreaFraction + GetHbvAquifer()->GetTotalHbvAreaFraction() / 100.0;
      elementFound = true;
    }
  if (GetKwaAquifer())
    {
      sumAreaFraction = sumAreaFraction + GetKwaAquifer()->GetTotalKwaAreaFraction() / 100.0;
      elementFound = true;
    }
  if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
      cout << " sumAreaFraction GetLandArea " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
      exit(1);
    }
  if (elementFound)
    {
      return sumAreaFraction * GetArea();
    }
  else
    {
      return missingData;
    }
}

// Algorithm to be performed in case: no input to landscape element from upstream elements
//void DistributedElement::WaterBalance(int timeStep) const
// Algorithm to be performed in case: input to landscape element from upstream elements
void DistributedElement::WaterBalance(int timeStep, int initialTimeSteps, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge) const
{
  HbvAquifer *lastHbvAquifer;
  KwaAquifer *lastKwaAquifer;
  if (GetLake()) {
    GetLake()->WaterBalance(timeStep, upLandAccumulatedLowerDischarge+upLandAccumulatedUpperDischarge);
  }
  if (GetGlacier()) {
    GetGlacier()->WaterBalance(timeStep, initialTimeSteps, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
  }
  if (GetHbvAquifer()) {
    lastHbvAquifer=GetHbvAquifer();
    while (lastHbvAquifer) {
      // Algorithm to be performed in case: no input to landscape element from upstream elements
      //      lastHbvAquifer->WaterBalance(timeStep);
      // Algorithm to be performed in case: input to landscape element from upstream elements
      lastHbvAquifer->WaterBalance(timeStep, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
      lastHbvAquifer=lastHbvAquifer->GetNextHbvAquifer();
    }
  }
  if (GetKwaAquifer()) {
    lastKwaAquifer=GetKwaAquifer();
    while (lastKwaAquifer) {
      lastKwaAquifer->WaterBalance(timeStep, upLandAccumulatedLowerDischarge, upLandAccumulatedUpperDischarge);
      lastKwaAquifer=lastKwaAquifer->GetNextKwaAquifer();
    }
  }
}

double DistributedElement::GetPrecipitation() const
{
  double precipitation=0.0;
  if (GetLake()) precipitation=GetLake()->GetPrecipitation();
  if (GetGlacier()) precipitation=precipitation+GetGlacier()->GetPrecipitation();
  if (GetHbvAquifer()) precipitation=precipitation+GetHbvAquifer()->GetPrecipitation();
  if (GetKwaAquifer()) precipitation=precipitation+GetKwaAquifer()->GetPrecipitation();
  return precipitation;
}

double DistributedElement::GetTemperature() const
{
  double temperature=0.0;
  if (GetLake()) temperature=GetLake()->GetTemperature();
  if (GetGlacier()) temperature=temperature+GetGlacier()->GetTemperature();
  if (GetHbvAquifer()) temperature=temperature+GetHbvAquifer()->GetTemperature();
  if (GetKwaAquifer()) temperature=temperature+GetKwaAquifer()->GetTemperature();
  return temperature;
}

double DistributedElement::GetSnowStore() const
{
  double snowStore=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    snowStore=GetGlacier()->GetSnowStore();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    snowStore=snowStore+GetHbvAquifer()->GetSnowStore();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetKwaAquifer()) {
    snowStore=snowStore+GetKwaAquifer()->GetSnowStore();
    sumAreaFraction=sumAreaFraction+GetKwaAquifer()->GetTotalKwaAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetSnowStore " << elementFound << "  " << snowStore;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetSnowStore " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return snowStore/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetSnowCoverFraction() const
{
  double snowCoverFraction=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    snowCoverFraction=GetGlacier()->GetSnowCoverFraction();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    snowCoverFraction=snowCoverFraction+GetHbvAquifer()->GetSnowCoverFraction();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetKwaAquifer()) {
    snowCoverFraction=snowCoverFraction+GetKwaAquifer()->GetSnowCoverFraction();
    sumAreaFraction=sumAreaFraction+GetKwaAquifer()->GetTotalKwaAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetSnowCoverFraction " << elementFound << "  " << snowCoverFraction;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetSnowCoverFraction " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return snowCoverFraction/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetMeltWater() const
{
  double meltWater=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    meltWater=GetGlacier()->GetMeltWater();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    meltWater=meltWater+GetHbvAquifer()->GetMeltWater();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetKwaAquifer()) {
    meltWater=meltWater+GetKwaAquifer()->GetMeltWater();
    sumAreaFraction=sumAreaFraction+GetKwaAquifer()->GetTotalKwaAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetMeltWater " << elementFound << "  " << meltWater;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetMeltWater " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return meltWater/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetGlacierMassBalance() const
{
    double glacierMassBalance = 0.0;
    double sumAreaFraction = 0.0;
    bool elementFound = false;
    //  cout << "start GetGlacierMassBalance\n";
    if (GetGlacier() && GetGlacier()->GetAreaFraction() > 0.0)
    {
        glacierMassBalance = GetGlacier()->GetGlacierMassBalance();
        //    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
        sumAreaFraction = sumAreaFraction + GetGlacier()->GetAreaFraction() / 100.0;
        elementFound = true;
    }
    //  cout << " sumAreaFraction GetGlacierMassBalance " << elementFound << "  " << glacierMassBalance;
    //  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
    if (elementFound && (sumAreaFraction <= 0.0 || sumAreaFraction > 1.0 + epsilon))
    {
        cout << " sumAreaFraction GetGlacierMassBalance " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
        exit(1);
    }
    //  cout << "  end GetGlacierMassBalance\n";
    //  if (elementFound) return 1.0;
    //  if (elementFound) return glacierMassBalance/sumAreaFraction;
    if (elementFound)
    {
        return glacierMassBalance/sumAreaFraction;
	//        return glacierMassBalance;
    }
    //  if (elementFound) return glacierMassBalance;
    else
    {
        return missingData;
    }
}

double DistributedElement::GetKiWaSoilMoisture(double lengthFraction) const
{
  double soilMoisture=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    soilMoisture=GetGlacier()->GetKiWaSoilMoisture(lengthFraction);
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetKwaAquifer()) {
    soilMoisture=soilMoisture+GetKwaAquifer()->GetSoilMoisture(lengthFraction);
    sumAreaFraction=sumAreaFraction+GetKwaAquifer()->GetTotalKwaAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetKiWaSoilMoisture " << elementFound << "  " << soilMoisture;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetKiWaSoilMoisture " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return soilMoisture/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetKiWaGroundWaterDepth(double lengthFraction) const
{
  double groundWaterDepth=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    groundWaterDepth=GetGlacier()->GetKiWaGroundWaterDepth(lengthFraction);
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetKwaAquifer()) {
    groundWaterDepth=groundWaterDepth+GetKwaAquifer()->GetGroundWaterDepth(lengthFraction);
    sumAreaFraction=sumAreaFraction+GetKwaAquifer()->GetTotalKwaAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetKiWaGroundWaterDepth " << elementFound << "  " << groundWaterDepth;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetKiWaGroundWaterDepth " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return groundWaterDepth/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetHbvSoilMoisture() const
{
  double soilMoisture=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    soilMoisture=GetGlacier()->GetHbvSoilMoisture();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    soilMoisture=soilMoisture+GetHbvAquifer()->GetSoilMoisture();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetSoilMoisture " << elementFound << "  " << soilMoisture;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetSoilMoisture " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return soilMoisture/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetHbvSoilMoistureDeficit() const
{
  double soilMoistureDeficit=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    soilMoistureDeficit=GetGlacier()->GetHbvSoilMoistureDeficit();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    soilMoistureDeficit=soilMoistureDeficit+GetHbvAquifer()->GetSoilMoistureDeficit();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetSoilMoistureDeficit " << elementFound << "  " << soilMoistureDeficit;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetSoilMoistureDeficit " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return soilMoistureDeficit/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetHbvPercSoilUpper() const
{
  double percSoilUpper=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    percSoilUpper=GetGlacier()->GetHbvPercSoilUpper();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    percSoilUpper=percSoilUpper+GetHbvAquifer()->GetPercSoilUpper();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetPercSoilUpper " << elementFound << "  " << percSoilUpper;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetPercSoilUpper " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return percSoilUpper/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetHbvUpperZone() const
{
  double upperZone=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    upperZone=GetGlacier()->GetHbvUpperZone();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    upperZone=upperZone+GetHbvAquifer()->GetUpperZone();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetHbvUpperZone " << elementFound << "  " << upperZone;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetHbvUpperZone " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return upperZone/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetHbvLowerZone() const
{
  double lowerZone=0.0;
  double sumAreaFraction=0.0;
  bool elementFound=false;
  if (GetGlacier()) {
    lowerZone=GetGlacier()->GetHbvLowerZone();
    sumAreaFraction=sumAreaFraction+GetGlacier()->GetAreaFraction()/100.0;
    elementFound=true;
  }
  if (GetHbvAquifer()) {
    lowerZone=lowerZone+GetHbvAquifer()->GetLowerZone();
    sumAreaFraction=sumAreaFraction+GetHbvAquifer()->GetTotalHbvAreaFraction()/100.0;
    elementFound=true;
  }
//  cout << " sumAreaFraction GetHbvLowerZone " << elementFound << "  " << lowerZone;
//  cout << "  " << GetLandIndex() << "  " << sumAreaFraction << endl;
  if (elementFound && (sumAreaFraction<=0.0 || sumAreaFraction>1.0)) {
    cout << " sumAreaFraction GetHbvLowerZone " << GetLandIndex() << "  " << sumAreaFraction << endl << endl;
    exit(1);
  }
  if (elementFound) return lowerZone/sumAreaFraction;
  else return missingData;
}

double DistributedElement::GetInterceptionLoss() const
{
  double interceptionLoss=0.0;
  if (GetHbvAquifer()) interceptionLoss=GetHbvAquifer()->GetInterceptionLoss();
  if (GetKwaAquifer()) interceptionLoss=interceptionLoss+GetKwaAquifer()->GetInterceptionLoss();
//  cout << " interceptionLoss " << interceptionLoss << endl;
  return interceptionLoss;
}

double DistributedElement::GetTranspSoilEvap() const
{
  double transpSoilEvap=0.0;
  if (GetHbvAquifer()) transpSoilEvap=GetHbvAquifer()->GetTranspSoilEvap();
  if (GetKwaAquifer()) transpSoilEvap=transpSoilEvap+GetKwaAquifer()->GetTranspSoilEvap();
  return transpSoilEvap;
}

double DistributedElement::GetLakeEvap() const
{
  double lakeEvap=0.0;
  if (GetLake()) lakeEvap=GetLake()->GetLakeEvap();
  return lakeEvap;
}

double DistributedElement::GetRunoff() const
{
  double runoff=0.0;
  if (GetLake()) runoff=GetLake()->GetRunoff();
  if (GetGlacier()) runoff=runoff+GetGlacier()->GetRunoff();
  if (GetHbvAquifer()) runoff=runoff+GetHbvAquifer()->GetRunoff();
  if (GetKwaAquifer()) runoff=runoff+GetKwaAquifer()->GetRunoff();
  return runoff;
}

double DistributedElement::GetLowerRunoff() const
{
  double lowerRunoff=0.0;
  if (GetLake()) lowerRunoff=GetLake()->GetRunoff();
  if (GetGlacier()) lowerRunoff=lowerRunoff+GetGlacier()->GetLowerRunoff();
  if (GetHbvAquifer()) lowerRunoff=lowerRunoff+GetHbvAquifer()->GetLowerRunoff();
  if (GetKwaAquifer()) lowerRunoff=lowerRunoff+GetKwaAquifer()->GetLowerRunoff();
  return lowerRunoff;
}

double DistributedElement::GetUpperRunoff() const
{
  double upperRunoff=0.0;
  if (GetGlacier()) upperRunoff=upperRunoff+GetGlacier()->GetUpperRunoff();
  if (GetHbvAquifer()) upperRunoff=upperRunoff+GetHbvAquifer()->GetUpperRunoff();
  if (GetKwaAquifer()) upperRunoff=upperRunoff+GetKwaAquifer()->GetUpperRunoff();
  return upperRunoff;
}

double DistributedElement::GetDischarge() const
{
  return GetRunoff()*area/commonPar->GetSECONDS_TIMESTEP();              /*  runoff (m) -> discharge (m3/s)  */
}

double DistributedElement::GetLowerDischarge() const
{
  return GetLowerRunoff()*area/commonPar->GetSECONDS_TIMESTEP();         /*  runoff (m) -> discharge (m3/s)  */
}

double DistributedElement::GetUpperDischarge() const
{
  return GetUpperRunoff()*area/commonPar->GetSECONDS_TIMESTEP();         /*  runoff (m) -> discharge (m3/s)  */
}

void DistributedElement::SetTotalReservoirStorage(DateTimeInformation * const DateTimeStore, int index, bool * firstTotal, bool inputDataFound)
{
  int initialTimeSteps;
  double totalStorage;
  initialTimeSteps = DateTimeStore->GetInitialTimeSteps();
  // Initial value
  if (*firstTotal) {
    if (index > 0) totalStorage = totalReservoir->GetTotalReservoirStorage(index-1);
    else totalStorage = totalReservoir->GetInitialTotalReservoirStorage();
  }
  else {
    totalStorage = totalReservoir->GetTotalReservoirStorage(index);
  }
  //  if (inputDataFound && index >= initialTimeSteps) {
  if (inputDataFound) {
    if (GetLake()) totalStorage = totalStorage + GetLake()->GetLakeStorageChange()*area;
  }
  totalReservoir->SetTotalReservoirStorage(index, totalStorage);
}


// class WaterCourse
class WaterCourse
{
public:
  WaterCourse();
  ~WaterCourse(); 
  void SetWaterCourseIndex(int value) { waterCourseIndex = value; }
  int GetWaterCourseIndex() const { return waterCourseIndex; }
  void SetIdentifier(int value) { identifier = value; }
  int GetIdentifier() const { return identifier; }
  void SetNumUpStream(int value);
  int GetNumUpStream() const { return numUpStream; }
  void SetUpStream(int k, WaterCourse *theWaterCourse) { UpStream[k] = theWaterCourse; }
  WaterCourse *GetUpStream(int k) const { return UpStream[k]; }
  void SetNumLandScape(int value) { numLandScape = value; }
  int GetNumLandScape() const { return numLandScape; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }
  void SetLakeNumber(int value) { lakeNumber = value; }
  int GetLakeNumber() const { return lakeNumber; }
  void SetSelectedWaterCourseTimeSeriesElements(SelectedWaterCourseTimeSeriesElements *object) { selectedWaterCourseElements = object; }
  SelectedWaterCourseTimeSeriesElements *GetSelectedWaterCourseTimeSeriesElements() const { return selectedWaterCourseElements; }
  void SetCorrection(double value) { correction = value; }
  double GetCorrection() const { return correction; }
  void AllocateAccumulatedDischarge(int numberTimeSteps); 
  void SetAccumulatedDischarge(int index, double value) { accumulatedDischarge[index] = value; }
  double GetAccumulatedDischarge(int index) const { return accumulatedDischarge[index]; }
  void AllocateAccumulatedInFlow(int numberTimeSteps); 
  void SetAccumulatedInFlow(int index, double value) { accumulatedInFlow[index] = value; }
  double GetAccumulatedInFlow(int index) const { return accumulatedInFlow[index]; }
  void AllocateAccumulatedWaterBalance(int numberTimeSteps); 
  void AllocateWaterBalance(int numberTimeSteps); 
  void SetAccumulatedPrecipitation(int index, double value) { accumulatedPrecipitation[index] = value; }
  double GetAccumulatedPrecipitation(int index) const { return accumulatedPrecipitation[index]; }
  void SetAccumulatedTemperature(int index, double value) { accumulatedTemperature[index] = value; }
  double GetAccumulatedTemperature(int index) const { return accumulatedTemperature[index]; }
  void SetAccumulatedSnowStore(int index, double value) { accumulatedSnowStore[index] = value; }
  double GetAccumulatedSnowStore(int index) const { return accumulatedSnowStore[index]; }
  void SetAccumulatedGlacierMassBalance(int index, double value) { accumulatedGlacierMassBalance[index] = value; }
  double GetAccumulatedGlacierMassBalance(int index) const { return accumulatedGlacierMassBalance[index]; }
  void SetAccumulatedEvapotranspiration(int index, double value) { accumulatedEvapotranspiration[index] = value; }
  double GetAccumulatedEvapotranspiration(int index) const { return accumulatedEvapotranspiration[index]; }
  void SetAccumulatedRunoff(int index, double value) { accumulatedRunoff[index] = value; }
  double GetAccumulatedRunoff(int index) const { return accumulatedRunoff[index]; }
  void SetAccumulatedHbvSoilMoisture(int index, double value) { accumulatedHbvSoilMoisture[index] = value; }
  double GetAccumulatedHbvSoilMoisture(int index) const { return accumulatedHbvSoilMoisture[index]; }
  void SetAccumulatedHbvSoilMoistureDeficit(int index, double value) { accumulatedHbvSoilMoistureDeficit[index] = value; }
  double GetAccumulatedHbvSoilMoistureDeficit(int index) const { return accumulatedHbvSoilMoistureDeficit[index]; }
  void SetAccumulatedHbvPercSoilUpper(int index, double value) { accumulatedHbvPercSoilUpper[index] = value; }
  double GetAccumulatedHbvPercSoilUpper(int index) const { return accumulatedHbvPercSoilUpper[index]; }
  void SetAccumulatedHbvUpperZone(int index, double value) { accumulatedHbvUpperZone[index] = value; }
  double GetAccumulatedHbvUpperZone(int index) const { return accumulatedHbvUpperZone[index]; }
  void SetAccumulatedHbvLowerZone(int index, double value) { accumulatedHbvLowerZone[index] = value; }
  double GetAccumulatedHbvLowerZone(int index) const { return accumulatedHbvLowerZone[index]; }
  void SetAccumulatedSum(int index, double value) { accumulatedSum[index] = value; }
  double GetAccumulatedSum(int index) const { return accumulatedSum[index]; }
  void SetAccumulatedSumSnow(int index, double value) { accumulatedSumSnow[index] = value; }
  double GetAccumulatedSumSnow(int index) const { return accumulatedSumSnow[index]; }
  void SetAccumulatedSumGlacier(int index, double value) { accumulatedSumGlacier[index] = value; }
  double GetAccumulatedSumGlacier(int index) const { return accumulatedSumGlacier[index]; }
  void SetAccumulatedSumHbv(int index, double value) { accumulatedSumHbv[index] = value; }
  double GetAccumulatedSumHbv(int index) const { return accumulatedSumHbv[index]; }
  void SetWaterCoursePrecipitation(int index, double value) { waterCoursePrecipitation[index] = value; }
  double GetWaterCoursePrecipitation(int index) const { return waterCoursePrecipitation[index]; }
  void SetWaterCourseTemperature(int index, double value) { waterCourseTemperature[index] = value; }
  double GetWaterCourseTemperature(int index) const { return waterCourseTemperature[index]; }
  void SetWaterCourseSnowStore(int index, double value) { waterCourseSnowStore[index] = value; }
  double GetWaterCourseSnowStore(int index) const { return waterCourseSnowStore[index]; }
  void SetWaterCourseGlacierMassBalance(int index, double value) { waterCourseGlacierMassBalance[index] = value; }
  double GetWaterCourseGlacierMassBalance(int index) const { return waterCourseGlacierMassBalance[index]; }
  void SetWaterCourseEvapotranspiration(int index, double value) { waterCourseEvapotranspiration[index] = value; }
  double GetWaterCourseEvapotranspiration(int index) const { return waterCourseEvapotranspiration[index]; }
  void SetWaterCourseRunoff(int index, double value) { waterCourseRunoff[index] = value; }
  double GetWaterCourseRunoff(int index) const { return waterCourseRunoff[index]; }
  void SetWaterCourseHbvSoilMoisture(int index, double value) { waterCourseHbvSoilMoisture[index] = value; }
  double GetWaterCourseHbvSoilMoisture(int index) const { return waterCourseHbvSoilMoisture[index]; }
  void SetWaterCourseHbvSoilMoistureDeficit(int index, double value) 
    { waterCourseHbvSoilMoistureDeficit[index] = value; }
  double GetWaterCourseHbvSoilMoistureDeficit(int index) const 
    { return waterCourseHbvSoilMoistureDeficit[index]; }
  void SetWaterCourseHbvPercSoilUpper(int index, double value) { waterCourseHbvPercSoilUpper[index] = value; }
  double GetWaterCourseHbvPercSoilUpper(int index) const { return waterCourseHbvPercSoilUpper[index]; }
  void SetWaterCourseHbvUpperZone(int index, double value) { waterCourseHbvUpperZone[index] = value; }
  double GetWaterCourseHbvUpperZone(int index) const { return waterCourseHbvUpperZone[index]; }
  void SetWaterCourseHbvLowerZone(int index, double value) { waterCourseHbvLowerZone[index] = value; }
  double GetWaterCourseHbvLowerZone(int index) const { return waterCourseHbvLowerZone[index]; }
  void ObsDataInput(DateTime firstTime, DateTime lastTime, int numberTimeSteps, int secondsPerTimeStep);
  double GetObsData(int index) const { return observedData[index]; }

private:
  int waterCourseIndex;
  int identifier;
  int lakeNumber;
  int numUpStream;
  int numLandScape;
  double correction;
  double *accumulatedDischarge;
  double *accumulatedInFlow;
  double *accumulatedPrecipitation;
  double *accumulatedTemperature;
  double *accumulatedSnowStore;
  double *accumulatedGlacierMassBalance;
  double *accumulatedEvapotranspiration;
  double *accumulatedRunoff;
  double *accumulatedHbvSoilMoisture;
  double *accumulatedHbvSoilMoistureDeficit;
  double *accumulatedHbvPercSoilUpper;
  double *accumulatedHbvUpperZone;
  double *accumulatedHbvLowerZone;
  double *accumulatedSum;
  double *accumulatedSumSnow;
  double *accumulatedSumGlacier;
  double *accumulatedSumHbv;
  double *waterCoursePrecipitation;
  double *waterCourseTemperature;
  double *waterCourseSnowStore;           
  double *waterCourseGlacierMassBalance;
  double *waterCourseEvapotranspiration;  
  double *waterCourseRunoff;              
  double *waterCourseHbvSoilMoisture;     
  double *waterCourseHbvSoilMoistureDeficit;     
  double *waterCourseHbvPercSoilUpper;     
  double *waterCourseHbvUpperZone;        
  double *waterCourseHbvLowerZone;        
  double *observedData;
  WaterCourse **UpStream;
  DistributedElement *landScapeElement;
  SelectedWaterCourseTimeSeriesElements *selectedWaterCourseElements;
};

WaterCourse::WaterCourse():
  lakeNumber(-9999),
  numUpStream(0),
  numLandScape(0),
  correction(0.0)
{
  UpStream = 0;
  SetLandScapeElement(0); 
  SetSelectedWaterCourseTimeSeriesElements(0);
}

WaterCourse::~WaterCourse()
{ 
}

void WaterCourse::SetNumUpStream(int value) 
{
  numUpStream = value;
  WaterCourse **upStr = new WaterCourse * [value];
  UpStream = upStr;
}

void WaterCourse::AllocateAccumulatedDischarge(int numberTimeSteps) 
{
  accumulatedDischarge = new double [numberTimeSteps];
}

void WaterCourse::AllocateAccumulatedInFlow(int numberTimeSteps) 
{
  accumulatedInFlow = new double [numberTimeSteps];
}

void WaterCourse::AllocateAccumulatedWaterBalance(int numberTimeSteps) 
{
  accumulatedPrecipitation = new double [numberTimeSteps];
  accumulatedTemperature = new double [numberTimeSteps];
  accumulatedSnowStore = new double [numberTimeSteps];
  accumulatedGlacierMassBalance = new double[numberTimeSteps];
  accumulatedEvapotranspiration = new double [numberTimeSteps];
  accumulatedRunoff = new double [numberTimeSteps];
  accumulatedHbvSoilMoisture = new double [numberTimeSteps];
  accumulatedHbvSoilMoistureDeficit = new double [numberTimeSteps];
  accumulatedHbvPercSoilUpper = new double [numberTimeSteps];
  accumulatedHbvUpperZone = new double [numberTimeSteps];
  accumulatedHbvLowerZone = new double [numberTimeSteps];
  accumulatedSum = new double [numberTimeSteps];
  accumulatedSumSnow = new double [numberTimeSteps];
  accumulatedSumGlacier = new double[numberTimeSteps];
  accumulatedSumHbv = new double [numberTimeSteps];
}

void WaterCourse::AllocateWaterBalance(int numberTimeSteps) 
{
  waterCoursePrecipitation = new double [numberTimeSteps];
  waterCourseTemperature = new double [numberTimeSteps];
  waterCourseSnowStore = new double [numberTimeSteps];
  waterCourseGlacierMassBalance = new double[numberTimeSteps];
  waterCourseEvapotranspiration = new double [numberTimeSteps];
  waterCourseRunoff = new double [numberTimeSteps];
  waterCourseHbvSoilMoisture = new double [numberTimeSteps];
  waterCourseHbvSoilMoistureDeficit = new double [numberTimeSteps];
  waterCourseHbvPercSoilUpper = new double [numberTimeSteps];
  waterCourseHbvUpperZone = new double [numberTimeSteps];
  waterCourseHbvLowerZone = new double [numberTimeSteps];
}

void WaterCourse::ObsDataInput(DateTime firstTime, DateTime lastTime, int numberTimeSteps, int secondsPerTimeStep)
{
  char c, buffer[256];
  int i, waterCourseId, year, month, day, hour, min;
  double value;
  bool catchmentFound=false;
  DateTime * datetime = new DateTime [numberTimeSteps];
  for (i=0; i<numberTimeSteps; i++) {
    datetime[i] = firstTime + i*secondsPerTimeStep;
  }
  if (datetime[0] != firstTime || datetime[numberTimeSteps-1] != lastTime) {
    cout << endl << " DateTime error during initialisation of observed data array: " <<
      datetime[0] << "  " << "  " << firstTime  << "  " << datetime[numberTimeSteps-1]  << "  " << lastTime << endl << endl;
    exit(1);
  }
  observedData = new double [numberTimeSteps];
  for (i=0; i<numberTimeSteps; i++) observedData[i] = missingData;
  ifstream fileIn("obs_streamflow.txt");
  if (fileIn == NULL) {
    cout << endl << " File obs_streamflow.txt not found" << endl << endl;
    //    exit (1);
  }
  while (fileIn.getline(buffer, 256) && !catchmentFound) {
    if (buffer[0]=='#') {
      sscanf(buffer,"%c %c %c %c %c %d ",&c,&c,&c,&c,&c,&waterCourseId);
      if (waterCourseId==GetIdentifier()) {
        i=0;
        while (fileIn.peek() != '#' && fileIn.peek() != EOF) {
          fileIn >> buffer >> value;
          sscanf(buffer,"%4d %2d %2d %c %2d %2d ", &year, &month, &day, &c, &hour, &min);
          while ((datetime[i].getYear()<year ||
                  (datetime[i].getYear()==year && datetime[i].getMonth()<month) ||
                  (datetime[i].getYear()==year && datetime[i].getMonth()==month && datetime[i].getDay()<day) ||
                  (datetime[i].getYear()==year && datetime[i].getMonth()==month && datetime[i].getDay()==day && 
                   datetime[i].getHour()<hour) ||
                  (datetime[i].getYear()==year && datetime[i].getMonth()==month && datetime[i].getDay()==day && 
                   datetime[i].getHour()==hour && datetime[i].getMinute()<min)) &&
                 i<numberTimeSteps-1) {
            i++;
            //      cout << " i = " << i << endl;
          }
          //    cout << datetime[i].getYear() << " " << datetime[i].getMonth() << " " << datetime[i].getDay() << " " << 
          //      datetime[i].getHour() << " " << datetime[i].getMinute() << endl;
          if (datetime[i].getYear()==year && datetime[i].getMonth()==month && datetime[i].getDay()==day && 
              datetime[i].getHour()==hour && datetime[i].getMinute()==min) {
            //      cout << " data found  " << i << " " << year << " " << month << " " << day << " " << hour << " " << min << " " << endl;
            if (value <= -9.0) 
              observedData[i] = missingData;
            else if (observedData[i] > -9.0 && observedData[i] < 0.0) 
              observedData[i] = 0.0;
            else
              observedData[i] = value;
            //    cout << waterCourseId << " " << i << " " << year << " " << month << " ";
            //    cout << day << " " << hour << " " << min << " " << " " << value << endl;
          }
          fileIn.ignore(256,'\n');
        }
      } else {
        while (fileIn.peek() != '#' && fileIn.peek() != EOF) {
          //      fileIn.getline(buffer, 256);
          fileIn.ignore(256,'\n');
        }
      }
    }
  }
  fileIn.close();
}  
