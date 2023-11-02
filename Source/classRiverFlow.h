/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments,          *
 *  river networks, lakes and landscape elements.                                       *
 *  Preprocessing for traversing and defining a hierarchy of watercourse elements.      *
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
#define ELEMENT(a,b) (((a)*nCols)+(b))

class DistributedElement
{
 public:
  DistributedElement();
  ~DistributedElement(); 
  void SetGeoIndex(int value) { geoIndex = value; }
  int GetGeoIndex() const { return geoIndex; }
  void SetLandIndex(int value) { landIndex = value; }
  int GetLandIndex() const { return landIndex; }
  void SetFlowIndex(int value) { flowIndex = value; }
  int GetFlowIndex() const { return flowIndex; }
  void SetFlowLandValue(int value) { flowDirectionLand = value; }
  int GetFlowLandValue() const { return flowDirectionLand; }
  void SetFlowRiverValue(int value) { flowDirectionRiver = value; }
  int GetFlowRiverValue() const { return flowDirectionRiver; }
  void SetLakeNumber(int value) { lakeNumber = value; }
  int GetLakeNumber() const { return lakeNumber; }
  void SetNumUpLand(int value) { numUpLand = value; }
  int GetNumUpLand() const { return numUpLand; }
  void SetWaterCourseElementIdentification(int value) { waterCourseElementIdentification = value; }
  int GetWaterCourseElementIdentification() const { return waterCourseElementIdentification; }
  void SetUpLandFlow(int k, DistributedElement *theElement) {upLandFlow[k] = theElement; }
  DistributedElement *GetUpLandFlow(int k) const { return upLandFlow[k]; }

private:
  int geoIndex;
  int waterCourseElementIdentification;
  int landIndex;
  int flowIndex;
  int flowDirectionLand;
  int flowDirectionRiver;
  int lakeNumber;
  int numUpLand;
  double area;
  DistributedElement *upLandFlow[7];
  DistributedElement *nextElement;
};

DistributedElement::DistributedElement()
{
  SetNumUpLand(0);
  for (int i=0; i<7; i++) 
    upLandFlow[i]=0;
}

DistributedElement::~DistributedElement()
{ 
}

class WaterCourseElement
{
 public:
  WaterCourseElement();
  ~WaterCourseElement();
  void SetGeoIndex(int value) { geoIndex = value; }
  int GetGeoIndex() const { return geoIndex; }
  void SetFlowIndex(int value) { flowIndex = value; }
  int GetFlowIndex() const { return flowIndex; }
  void SetNumWaterCourse(int value) { numWaterCourse = value; }
  int GetNumWaterCourse() const { return numWaterCourse; }
  void SetNumUpStream(int value) { numUpStream = value; }
  int GetNumUpStream() const { return numUpStream; }
  void SetUpStreamElement(int k, WaterCourseElement *theElement) {upStreamElement[k] = theElement; }
  WaterCourseElement *GetUpStream(int k) const { return upStreamElement[k]; }
  void SetLandScapeElement(DistributedElement *theElement) { landScapeElement = theElement; }
  DistributedElement *GetLandScapeElement() const { return landScapeElement; }
  void SetLakeElementList(DistributedElement *theElement) { lakeElementList = theElement; }
  DistributedElement *GetLakeElementList() const { return lakeElementList; }

private:
  int geoIndex;
  int flowIndex;
  int numWaterCourse;
  int numUpStream;
  WaterCourseElement *upStreamElement[7];
  DistributedElement *landScapeElement;
  DistributedElement *lakeElementList;
};

WaterCourseElement::WaterCourseElement()
{
  SetNumUpStream(0);
  SetLandScapeElement(0); 
  SetLakeElementList(0);
  for (int i=0; i<7; i++) 
    upStreamElement[i]=0;
}

WaterCourseElement::~WaterCourseElement()
{
}
