#======================================================================================================================!
#
#                    DassFlow Version 2.0
#
#======================================================================================================================!
#
#  Copyright University of Toulouse-INSA & CNRS (France)
#
#  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
#  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
#  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
#  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
#
#  DassFlow software includes few mostly independent "modules" with common architectures and structures:
#    - Shallow Module (Shallow Water Model, Finite Volume Method), i.e. the present code.
#    - 3D Module (Full Stokes Model, Finite Element Method, Mobile Gometries, ALE).
#  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
#
#  Many people have contributed to the DassFlow development from the initial version to the latest ones.
# 	Current main developer:
#               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
# 	with scientific and/or programming contributions of:
#               R. Madec   (Mathematics Institute of Toulouse IMT).
#               K. Larnier (Fluid Mechanics Institute of Toulouse IMFT).
#               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
#               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
#	and former other developers (M. Honnorat and J. Marin).
#
#  Scientific Contact : jerome.monnier@insa-toulouse.fr
#  Technical  Contact : frederic.couderc@math.univ-toulouse.fr
#
#  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
#  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
#  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
#
#  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
#  license, users are provided only with a limited warranty and the software's author, the holder of the economic
#  rights, and the successive licensors have only limited liability.
#
#  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
#  developing or reproducing the software by the user in light of its specific status of free software, that may
#  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
#  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
#  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
#  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
#  accept its terms.
#
#======================================================================================================================!

#!/usr/bin/env perl

use strict;
use warnings;

#-----------------------------------------------------------------------------------------------------------------------
# Replace function calc_cost_and_gradients in dassflow1d/__init__.py
open(my $in, "<dassflow1d/__init__.py"); # Ouvrir le fichier "file" en mode "lecture"

my @copy = <$in>;  #Copier le contenu du fichier dans "copy"

close $in; #Fermer l'original

open(my $out, ">dassflow1d/__init__.py"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

my $step = 0;
my $line = 0;

for ( @copy ) #Parcourt la copie du fichier ligne par ligne
{
    $line = $line + 1;
    if ( $line == 1 )
    {
        if ( m/from __future__(.*)/i )
        {
            print $out "$_";
            print $out "import numpy as np\n";
            next;
            
        } else {
            print $out "import numpy as np\n";
        }
    }
    if ( $step == 0)
    {
        if ( m/^(\s*)def calc_cost_and_gradients(.*)grad\, (.*)/i )
        {
            print $out "$1def calc_cost_and_gradients$2$3\n";
            $step = 1;
            next;
        }
        else
        {
          print $out "$_";
          next;
        }
    }
    elsif ( $step == 1 )
    {
        if ( m/^(\s*)cost = calc_cost_and_gradients(.*)grad\, (.*)/i )
        {
            print $out "$1cost, grad = calc_cost_and_gradients$2$3\n";
            next;
        }
        if ( m/^(\s*)grad \: float array.*/i )
        {
            next;
        }
        if ( m/^(\s*)cost \: float(.*)/i )
        {
            print $out "$_";
            print $out "$1grad \: float array$2\n";
            next;
        }
        if ( m/^(\s*)cost \= _dassflow1d(.*)/i )
        {
            print $out "$1grad = np.zeros(ctrl.x.size)\n$_";
            next;
        }
        if ( m/^(\s*)return cost(.*)/i )
        {
            print $out "$1return cost, grad$2\n";
            $step = 0;
            next;
        }
        else
        {
          print $out "$_";
          next;
        }
    }
}
        
# Add function calc_cost_and_gradients_scipy_minimize
print $out "\ndef calc_cost_and_gradients_scipy_minimize(x, mdl, ctrl, obs):\n";
print $out "\n    \"\"\"\"\n";
print $out "    cost, grad = calc_cost_and_gradients_scipy_minimize(x, mdl, ctrl, obs)\n";
print $out "\n\n";    
print $out "    Parameters\n";
print $out "    ----------\n";
print $out "    x : Current control vector values\n";
print $out "    mdl : Model\n";
print $out "    ctrl : Control\n";
print $out "    obs : Observations\n";
print $out "\n    Returns\n";
print $out "    -------\n";
print $out "    cost : float\n";
print $out "    grad : float array\n";
print $out "\n";
print $out "    \"\"\"\n\n";
print $out "    ctrl.x[:] = x[:]\n";
print $out "    grad = np.zeros(x.size)\n";
print $out "    cost = _dassflow1d.f90wrap_calc_cost_and_gradients(grad=grad, mdl=mdl._handle, \\\n";
print $out "        ctrl=ctrl._handle, obs=obs._handle)\n";
print $out "    mdl.__last_cost__ = cost\n";
print $out "    mdl.__last_grad__ = grad\n";
print $out "    return cost, grad\n";
print $out "\n";
close $out;

#-----------------------------------------------------------------------------------------------------------------------
# Add function get_item_slice in dassflow1d/m_control.py
open($in, "<dassflow1d/m_control.py"); # Ouvrir le fichier "file" en mode "lecture"

@copy = <$in>;  #Copier le contenu du fichier dans "copy"

close $in; #Fermer l'original

$step = 0;

open($out, ">dassflow1d/m_control.py"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

for ( @copy ) #Parcourt la copie du fichier ligne par ligne
{
    if ( $step == 0)
    {
        if (m/^(\s*)class Control\((.*)/i )
        {
            $step = 1;
        }
        print $out "$_";
        next;
    }
    else
    {
        if (m/^(\s*)def __str__\(self\):(.*)/i)
        {
            print $out "    def get_item_slice(self, index):\n";
            print $out "        if index < 0 or index > self.get_items_count():\n";
            print $out "            raise IndexError('index out of range')\n";
            print $out "        return slice(self.items[index].offset, self.items[index].offset+self.items[index].nx)\n";
            print $out "\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        else
        {
          if (m/^(\s*)class (.*)/i)
          {
              $step = 0;
          }
          print $out "$_";
          next;
        }
    }
}
close $out;

#-----------------------------------------------------------------------------------------------------------------------
# Add functions in dassflow1d/m_mesh.py
open($in, "<dassflow1d/m_mesh.py"); # Ouvrir le fichier "file" en mode "lecture"

@copy = <$in>;  #Copier le contenu du fichier dans "copy"

close $in; #Fermer l'original

$step = 0;

open($out, ">dassflow1d/m_mesh.py"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

for ( @copy ) #Parcourt la copie du fichier ligne par ligne
{
    if ( $step == 0)
    {
        if (m/^(\s*)class Mesh\((.*)/i )
        {
            $step = 1;
        }
        print $out "$_";
        next;
    }
    else
    {
        if (m/^(\s*)def __str__\(self\):(.*)/i)
        {
    
            print $out "    def get_segment_field(self, iseg, field):\n";
            print $out "    \n";
            print $out "        if field == \"curvilinear_abscissa\" or field == \"x\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            return [self.cs[i].x for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"bathy\" or field == \"b\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            return [self.cs[i].bathy for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"z0\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            return [self.cs[i].level_heights[0] for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"w0\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            return [self.cs[i].level_widths[0] for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"A0\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            w0 = [self.cs[i].level_widths[0] for i in range(first_cs-1, last_cs)]\n";
            print $out "            z0 = [self.cs[i].level_heights[0] for i in range(first_cs-1, last_cs)]\n";
            print $out "            bathy = [self.cs[i].bathy for i in range(first_cs-1, last_cs)]\n";
            print $out "            return [(z0[i]-bathy[i])*w0[i] for i in range(0, len(w0))]\n";
            print $out "        else:\n";
            print $out "          raise ValueError(\"Unknown segment field: \%s\" \% field)\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        else
        {
          if (m/^(\s*)class (.*)/i)
          {
              $step = 0;
          }
          print $out "$_";
          next;
        }
    }
}
close $out;

#-----------------------------------------------------------------------------------------------------------------------
# Add functions in dassflow1d/m_sw_mono.py
open($in, "<dassflow1d/m_sw_mono.py"); # Ouvrir le fichier "file" en mode "lecture"

@copy = <$in>;  #Copier le contenu du fichier dans "copy"

close $in; #Fermer l'original

$step = 0;

open($out, ">dassflow1d/m_sw_mono.py"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

for ( @copy ) #Parcourt la copie du fichier ligne par ligne
{
    if ( $step == 0)
    {
        if (m/^(\s*)class Unknowns\((.*)/i )
        {
            $step = 1;
        }
        if (m/^(\s*)class Model\((.*)/i )
        {
            $step = 2;
        }
        print $out "$_";
        next;
    }
    elsif ($step == 1)
    {
        if (m/^(\s*)def __str__\(self\):(.*)/i)
        {
    
            print $out "    def get_segment_field(self, mesh, iseg, field):\n";
            print $out "    \n";
            print $out "        if field == \"a\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            return [self.a[i] for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"q\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            return [self.q[i] for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"h\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            return [self.h[i] for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"z\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            return [self.h[i] + mesh.cs[i].bathy for i in range(first_cs-1, last_cs)]\n";
            print $out "        else:\n";
            print $out "          raise ValueError(\"Unknown segment field: \%s\" \% field)\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        else
        {
          if (m/^(\s*)class (.*)/i)
          {
              $step = 0;
          }
          print $out "$_";
          next;
        }
    }
    elsif ($step == 2)
    {
        if (m/^(\s*)def __str__\(self\):(.*)/i)
        {
            print $out "    def get_segment_results(self, iseg, iout, field):\n";
            print $out "    \n";
            print $out "        if field == \"q\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            return [self.res.q[i, iout] for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"h\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            return [self.res.h[i, iout] for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"z\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            return [self.res.h[i, iout] + self.msh.cs[i].bathy for i in range(first_cs-1, last_cs)]\n";
            print $out "        elif field == \"a\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            return [self.res.a[i, iout] for i in range(first_cs-1, last_cs)]\n";
            print $out "        else:\n";
            print $out "          raise ValueError(\"Unknown result field: \%s\" \% field)\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        else
        {
          if (m/^(\s*)class (.*)/i)
          {
              $step = 0;
          }
          print $out "$_";
          next;
        }
    }
}
close $out;


#-----------------------------------------------------------------------------------------------------------------------
open($out, ">env.sh"); # ouvrir le fichier en mode "ecriture"  (écrasement du fichier)

print $out "FILEPATH=`realpath \${BASH_SOURCE[0]}`\n";
print $out "DIR=`dirname \$FILEPATH`\n";
print $out "export PYTHONPATH=\$DIR\n";

close($out);

