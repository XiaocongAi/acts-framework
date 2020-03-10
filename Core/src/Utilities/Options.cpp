// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Utilities/Options.hpp"

#include <ostream>

namespace {
/// Helper function to print multiple elements in a container.
template <typename Iterator>
inline std::ostream&
printContainer(std::ostream& os,
               Iterator      begin,
               Iterator      end,
               const char*   separator)
{
  for (auto it = begin; it != end; ++it) {
    if (it != begin) { os << separator; }
    os << *it;
  }
  return os;
}
}  // namespace

std::ostream&
std::operator<<(std::ostream& os, const read_series& vec)
{
  return printContainer(os, vec.begin(), vec.end(), " ");
}

std::ostream&
std::operator<<(std::ostream& os, const read_range& vec)
{
  return printContainer(os, vec.begin(), vec.end(), " ");
}

std::ostream&
std::operator<<(std::ostream& os, const read_strings& vec)
{
  return printContainer(os, vec.begin(), vec.end(), " ");
}
