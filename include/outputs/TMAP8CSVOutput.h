#pragma once

#include "CSV.h"

class TMAP8CSVOutput : public CSV
{
public:
  static InputParameters validParams();
  TMAP8CSVOutput(const InputParameters & parameters);
};
