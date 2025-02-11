
/**
 * @file tpccalib/TpcSpaceChargeMatrixContainerv2.cc
 * @author Xudong Yu
 * @date January 2025
 * @brief Contains matrices needed for space charge trackbase reconstruction 1D (layer,radius,phi,z)
 */

#include "TpcSpaceChargeMatrixContainerv2.h"

//___________________________________________________________
TpcSpaceChargeMatrixContainerv2::TpcSpaceChargeMatrixContainerv2()
{
  // reset all matrix arrays
  TpcSpaceChargeMatrixContainerv2::Reset();
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::identify(std::ostream& out) const
{
  out << "TpcSpaceChargeMatrixContainerv2" << std::endl;
  out << "  bins: " << m_bins << std::endl;
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::get_grid_dimensions(int& bins) const
{
  bins = m_bins;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainerv2::get_grid_size() const
{
  return m_bins;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainerv2::get_cell_index(int i) const
{
  if (i < 0 || i >= m_bins)
  {
    return -1;
  }
  return i;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainerv2::get_entries(int cell_index) const
{
  // bound check
  if (!bound_check(cell_index))
  {
    return 0;
  }
  return m_entries[cell_index];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainerv2::get_lhs(int cell_index, int i, int j) const
{
  // bound check
  if (!bound_check(cell_index, i, j))
  {
    return 0;
  }
  return m_lhs[cell_index][get_flat_index(i, j)];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainerv2::get_rhs(int cell_index, int i) const
{
  // bound check
  if (!bound_check(cell_index, i))
  {
    return 0;
  }
  return m_rhs[cell_index][i];
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::Reset()
{
  // reset total number of bins
  const int totalbins = m_bins;

  // reset arrays
  m_entries = std::vector<int>(totalbins, 0);
  m_lhs = std::vector<matrix_t>(totalbins, {{}});
  m_rhs = std::vector<column_t>(totalbins, {{}});
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::set_grid_dimensions(int bins)
{
  m_bins = bins;
  Reset();
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_entries(int cell_index, int value)
{
  if (bound_check(cell_index))
  {
    m_entries[cell_index] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_lhs(int cell_index, int i, int j, float value)
{
  if (bound_check(cell_index, i, j))
  {
    m_lhs[cell_index][get_flat_index(i, j)] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_rhs(int cell_index, int i, float value)
{
  if (bound_check(cell_index, i))
  {
    m_rhs[cell_index][i] += value;
  }
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainerv2::add(const TpcSpaceChargeMatrixContainer& other)
{
  // check dimensions
  int bins = 0;
  other.get_grid_dimensions(bins);
  if ((m_bins != bins))
  {
    std::cout << "TpcSpaceChargeMatrixContainerv2::add - inconsistent grid sizes" << std::endl;
    return false;
  }

  // increment cell entries
  for (size_t cell_index = 0; cell_index < m_lhs.size(); ++cell_index)
  {
    add_to_entries(cell_index, other.get_entries(cell_index));
  }

  // increment left hand side matrices
  for (size_t cell_index = 0; cell_index < m_lhs.size(); ++cell_index)
  {
    for (int i = 0; i < m_ncoord; ++i)
    {
      for (int j = 0; j < m_ncoord; ++j)
      {
        add_to_lhs(cell_index, i, j, other.get_lhs(cell_index, i, j));
      }
    }
  }

  // increment right hand side matrices
  for (size_t cell_index = 0; cell_index < m_lhs.size(); ++cell_index)
  {
    for (int i = 0; i < m_ncoord; ++i)
    {
      add_to_rhs(cell_index, i, other.get_rhs(cell_index, i));
    }
  }

  return true;
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainerv2::bound_check(int cell_index) const
{
  if (cell_index < 0 || cell_index >= (int) m_rhs.size())
  {
    return false;
  }
  return true;
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainerv2::bound_check(int cell_index, int i) const
{
  if (cell_index < 0 || cell_index >= (int) m_rhs.size())
  {
    return false;
  }
  if (i < 0 || i >= m_ncoord)
  {
    return false;
  }
  return true;
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainerv2::bound_check(int cell_index, int i, int j) const
{
  if (cell_index < 0 || cell_index >= (int) m_lhs.size())
  {
    return false;
  }
  if (i < 0 || i >= m_ncoord)
  {
    return false;
  }
  if (j < 0 || j >= m_ncoord)
  {
    return false;
  }
  return true;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainerv2::get_flat_index(int i, int j) const
{
  return j + i * m_ncoord;
}
