/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments,          *
 *  river networks, lakes and landscape elements.                                       *
 *  Preprocessing.                                                                      *
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
const int numberLandSurfaceClasses=7;        // All possible land surface types, excluding lakes and glaciers
const int numberSoilClasses=7;               // All possible soil/subsurface types, excluding lakes and glaciers
const int finalYear=2100;
const int finalMonth=12;
const int finalDay=31;
const int finalHour=23;
const int finalMinute=59;
const double missingData=-9999.0;


//class ParametersCommon
class ParametersCommon
{
 public:
  ParametersCommon();
  ~ParametersCommon();
  void SetNUM_PREC_SERIES(int value) { NUM_PREC_SERIES = value; }
  int GetNUM_PREC_SERIES() const { return NUM_PREC_SERIES; }
  void SetNUM_TEMP_SERIES(int value) { NUM_TEMP_SERIES = value; }
  int GetNUM_TEMP_SERIES() const { return NUM_TEMP_SERIES; }
  
 private:
  int NUM_PREC_SERIES;
  int NUM_TEMP_SERIES;
};

ParametersCommon::ParametersCommon():
  NUM_PREC_SERIES(0),
  NUM_TEMP_SERIES(0)
{
}
     
ParametersCommon::~ParametersCommon()
{     
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


//class DistributedElement
class DistributedElement
{
 public:
  DistributedElement();
  ~DistributedElement(); 
  void SetGeoIndex(int value) { geoIndex = value; }
  int GetGeoIndex() const { return geoIndex; }
  void SetLandIndex(int value) { landIndex = value; }
  int GetLandIndex() const { return landIndex; }
  void SetFlowLandValue(int value) { flowDirectionLand = value; }
  int GetFlowLandValue() const { return flowDirectionLand; }
  void SetLakeNumber(int value) { lakeNumber = value; }
  int GetLakeNumber() const { return lakeNumber; }
  void SetWaterCourseValue(int value) { waterCourseValue = value; }
  int GetWaterCourseValue() const { return waterCourseValue; }
  void SetNumUpLand(int value) { numUpLand = value; }
  int GetNumUpLand() const { return numUpLand; }
  void SetArea(double value) { area = value; }
  double GetArea() const { return area; }
  void SetElevation(double value) { elevation = value; }
  double GetElevation() const { return elevation; }
  void SetSlopeLength(double value) { slopeLength = value; }
  double GetSlopeLength() const { return slopeLength; }
  void SetSlopeAngle(double value) { slopeAngle = value; }
  double GetSlopeAngle() const { return slopeAngle; }
  void SetAspect(double value) { aspect = value; }
  double GetAspect() const { return aspect; }
  void SetLakePercent(double value) { lakePercent = value; }
  double GetLakePercent() const { return lakePercent; }
  void SetForestPercent(double value) { forestPercent = value; }
  double GetForestPercent() const { return forestPercent; }
  void SetBogPercent(double value) { bogPercent = value; }
  double GetBogPercent() const { return bogPercent; }
  void SetGlacierPercent(double value) { glacierPercent = value; }
  double GetGlacierPercent() const { return glacierPercent; }
  void SetOpenLandPercent(double value) { openLandPercent = value; }
  double GetOpenLandPercent() const { return openLandPercent; }
  void SetAlpinePercent(double value) { alpinePercent = value; }
  double GetAlpinePercent() const { return alpinePercent; }
  void SetHeatherPercent(double value) { heatherPercent = value; }
  double GetHeatherPercent() const { return heatherPercent; }
  void SetBedrockPercent(double value) { bedrockPercent = value; }
  double GetBedrockPercent() const { return bedrockPercent; }
  void SetTreeLevel(double value) { treeLevel = value; }
  double GetTreeLevel() const { return treeLevel; }
  void SetXCoord(double value) { xCoord = value; }
  double GetXCoord() const { return xCoord; }
  void SetYCoord(double value) { yCoord = value; }
  double GetYCoord() const { return yCoord; }
  void SetUpLandFlow(int k, DistributedElement *theElement) {upLandFlow[k] = theElement; }
  DistributedElement *GetUpLandFlow(int k) const { return upLandFlow[k]; }
  void SetNextElement(DistributedElement *theElement) { nextElement = theElement; }
  DistributedElement *GetNextElement() const { return nextElement; }

private:
  int geoIndex;
  int landIndex;
  int flowDirectionLand;
  int waterCourseValue;
  int lakeNumber;
  int numUpLand;
  double area;
  double elevation;
  double slopeLength;
  double slopeAngle;
  double aspect;
  double lakePercent;
  double forestPercent;
  double bogPercent;
  double glacierPercent;
  double openLandPercent;
  double alpinePercent;
  double heatherPercent;
  double bedrockPercent;
  double treeLevel;
  double xCoord;
  double yCoord;
  DistributedElement *upLandFlow[7];
  DistributedElement *nextElement;
};

DistributedElement::DistributedElement():
  lakePercent(0.0),
  forestPercent(0.0),
  bogPercent(0.0),
  glacierPercent(0.0),
  openLandPercent(0.0),
  alpinePercent(0.0),
  heatherPercent(0.0),
  bedrockPercent(0.0),
  xCoord(0.0),
  yCoord(0.0)
{
  SetNextElement(0);
  SetNumUpLand(0);
  for (int i=0; i<7; i++) 
    upLandFlow[i]=0;
}

DistributedElement::~DistributedElement()
{ 
}


//class WaterCourse
class WaterCourse
{
public:
  WaterCourse();
  ~WaterCourse(); 
  void SetWaterCourseIndex(int value) { waterCourseIndex = value; }
  int GetWaterCourseIndex() const { return waterCourseIndex; }
  void SetIdentifier(int value) { identifier = value; }
  int GetIdentifier() const { return identifier; }
  void SetLakeNumber(int value) { lakeNumber = value; }
  int GetLakeNumber() const { return lakeNumber; }
  void SetNumLandScape(int value) { numLandScape = value; }
  int GetNumLandScape() const { return numLandScape; }
  void SetCorrection(double value) { correction = value; }
  double GetCorrection() const { return correction; }
  void SetNumUpStream(int value);
  int GetNumUpStream() const { return numUpStream; }
  void SetUpStream(int k, WaterCourse *theWaterCourse) { UpStream[k] = theWaterCourse; }
  WaterCourse *GetUpStream(int k) const { return UpStream[k]; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }

private:
  int waterCourseIndex, identifier;
  int lakeNumber;
  int numLandScape;
  int numUpStream;
  double correction;
  WaterCourse **UpStream;
  DistributedElement *landScapeElement;
};

WaterCourse::WaterCourse():
  numLandScape(0),
  numUpStream(0),
  correction(0.0)
{
  UpStream = 0;
  SetLandScapeElement(0); 
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
