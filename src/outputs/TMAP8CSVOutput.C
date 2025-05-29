#include "TMAP8CSVOutput.h"

registerMooseObject("TMAP8App", TMAP8CSVOutput);

InputParameters
TMAP8CSVOutput::validParams()
{
  return CSV::validParams();  // Don't add anything new
}

TMAP8CSVOutput::TMAP8CSVOutput(const InputParameters & parameters)
  : CSV(parameters)
{
}
