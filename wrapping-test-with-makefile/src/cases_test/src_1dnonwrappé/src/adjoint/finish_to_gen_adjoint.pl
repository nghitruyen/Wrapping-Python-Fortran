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
# Lecture de l'intégralité des fichiers du répertoire ./tap qui sont du type 'm_*_diff.f90' et 'm_*_back.f90'
#  - Le fichier est conservé si et seulement si il contient au moins une subroutine ou fonction différenciée
#  - Tous les fichiers ayant une subroutine différenciée sont agrégés à la liste  list
#  - Supprime les subroutines et les fonctions non différenciées
#-----------------------------------------------------------------------------------------------------------------------

my @list; #Initialisation de la variable list

foreach my $file ( <m_*_diff.f90> , <m_*_back.f90> ) # Parmi tous les fichiers de type 'm_*_diff.f90' et 'm_*_back.f90'
   {
   
   open(my $in, "<$file"); # Ouvrir le fichier "file" en mode "lecture"

   my @copy = <$in>;  #Copier le contenu du fichier dans "copy"

   close $in; #Fermer l'original

   open(my $out, ">$file"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

   my $test = 0; # Initialisation variable "test" à 0
   my $keep = 0; # Initialisation variable "keep" à 0
 
   for ( @copy ) #Parcourt la copie du fichier ligne par ligne
      { 

      if    ( m/^\s*SUBROUTINE (\w*)_(DIFF|BACK).*/i ) #Si la copie contient une subroutine differenciee 
         {
         if ( $1 !~ m/com_dof/i ) #Si la subroutine n'est pas com_dof
            {
            $keep = 1;              #On fixe keep à 1
            @list = (@list,$file);  #On ajoute le fichier à la liste
            }
         }

      elsif ( m/^\s*SUBROUTINE.*/i )  #Sinon si la copie contient une surboutine
         {
         $test = 1;     #On fixe test à 1
         print $out ""; 
         next;
         }
      
      #if    ( m/^\s*(real|integer).*FUNCTION\s*.*_(DIFF|BACK).*/i )  #Si la copie contient une fonction differenciee
     # if    ( m/^\s*[^'END'].*FUNCTION\s*.*_(DIFF|BACK).*/i )
      if    (( m/^\s*.*FUNCTION\s*.*_(DIFF|BACK).*/i )&&( m/^\s.[^'END']\s*.*/i))
         {
         $keep = 1;     #On fixe keep à 1
         }

      #elsif ( m/^\s*(real|integer).*FUNCTION\s*.*/i )  #Si la copie contient une fonction 
      #elsif ( m/^\s*[^'END'].*FUNCTION\s*.*/i )
      elsif (( m/^\s*.*FUNCTION\s*.*/i )&&( m/^\s.[^'END']\s*.*/i))
         {
         #print $_;
         $test = 1;     # On fixe test à 1
         print $out ""; 
         next;
         }

      if    ( m/^\s*END\s*SUBROUTINE.*_(DIFF|BACK).*/i ) #Si la copie contient une fin de subroutine differenciee
         {
         print $out $_; # On écrit la ligne dans le fichier out
         next;
         }

      elsif ( m/^\s*END\s*SUBROUTINE.*/i ) #Si la copie contient une fin de subroutine 
         {
         $test = 0;     #On fixe test à 0
         print $out ""; 
         next;
         }

      if    ( m/^\s*END\s*FUNCTION\s*.*_(DIFF|BACK).*/i ) #Si la copie  contient une fin de fonction differenciee
         { 
         print $out $_; # On écrit la ligne dans le fichier out
         next;
         }

      elsif ( m/^\s*END\s*FUNCTION\s*.*/i ) # Si la copie contient une fin de fonction
         {
         $test = 0;        # On fixe test à 0
         print $out "";    
         next;
         }

      if ( m/^\s*MODULE\s*([a-zA-Z0-9_]*)_(DIFF|BACK).*/i ) # Si la copie contient un module différencié
         {
         $test = 1;                                        # On fixe le test à 1
         print $out $_."\n  USE $1\n  USE M_TAP_VARS\n\n"; # On écrit la ligne de déf du module plus l'inclusion de mpi et tap_vars
         #print $out $_."\n  USE $1\n  USE M_MPI_$2\n  USE M_TAP_VARS\n\n";
         next;
         }

      if ( m/^\s*CONTAINS.*/i ) # Si la copie contient un contains
         {
         $test = 0; # On fixe test ) 0
         }

      if ( $test == 1 ) # Si la variable test est à 1
         {
         print $out ""; # On n'est écrit rien
         next;
         }

      print $out $_;  #On écrit la ligne dans le fichier out

      }

   close $out; # On ferme le fichier out

   if ( $keep == 0 ) #Si la variable keep est égal à 0 
      {
      unlink $file   #On supprime le fichier
      }

   }

#-----------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------

my (%saw,@list_sort)=();
undef %saw;
@list_sort = sort(grep(!$saw{$_}++, @list));

open(my $diffvars, ">diffvars.txt");   #Ouvre le fichier diffvars en mode écriture (écrasement des données)

foreach my $file ( <*_diff.f90> , <*_back.f90>,<*_back0.f90> )  # Parmi tous les fichiers de type '*_diff.f90' et '*_back.f90'
   {
   open(my $in, "<$file"); # Ouvrir le fichier "file" en mode "lecture"
   my @copy = <$in>;  #Copier le contenu du fichier dans "copy"
   close $in; #Fermer l'original
   open(my $out, ">$file"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

   for ( @copy ) #Parcourt la copie du fichier ligne par ligne
      {
      if ( m/^!.*Hint:\s*([\w%]*) .* \**([\w%]*)/i ) #Si l'occurence '!  Hint:' est trouvée
         {
         my $sub1 = $1; # Réccupération du tableau
         my $sub2 = $2; # Réccupération du nom du tableau

         if ( $sub1 !~ /ISIZE(.)OF.*/i ) #Si ce n'est pas un problème de taille de tableau
            {
            next; #On passe
            }

         # Corrected following line (fomerly evaluates to 1 has it is a match condition result, 
         # not extration
         # my $sub3 = ( $sub1 =~ /ISIZE(.)OF.*/i );
         my ($sub3) = ( $sub1 =~ /ISIZE(.)OF.*/i );
         print $diffvars "$file\t\t$sub1\t\t$sub2\n";

         for ( @copy ) 
            {
            if ( !( m/^!.*Hint:\s*([\w%]*) .* \**([\w%]*)/i )) #Si l'occurence '!  Hint:' n'est pas trouvée
               {
               s/$sub1\b/size($sub2,$sub3)/g;
               }
            }

         next;
         }
      # Exemple !  Hint: ISIZE1OFDrfdof_v should be the size of dimension 1 of array *dof%v
      #$1 : ISIZE1OFDrfdof_v
      #$2 : *dof%v




      if ( m/^\s*IMPLICIT NONE.*/i && $file !~ /m_.*/ ) # Recherche des IMPLICIT NONE dans les fichiers "m_..."
         {
         print $out "\n  USE M_TAP_VARS ! Added by Perl Script -> Need to be filled !!!\n\n".$_; # Remplacement par USE M_TAP_VARS
         next;
         }

      if ( m/^(\s*)USE (.*)_(DIFF|BACK)/i && $file !~ /m_.*/ )  # Recherche des USE ..._DIFF ou ..._BACK dans les fichiers différents de "m_..."
         {
         #print "$file : $_ ";
         print $out "$1USE $2 ! Replaced by Perl Script\n"; # Ajout de USE NOM module d'origine

         my $sub1 = $1;
         my $sub2 = $2;
         my $sub3 = $3;

         if ( "@list_sort" =~ /.*$sub2.*/i ) 
            {
            print $out $sub1 ; print $out "USE $sub2\_$sub3\n";
            }

         next;

         }

      if ( m/^(\s*)TYPE\((\w*)_(DIFF|BACK)\)(.*)/i ) 
         {
         print $out "$1TYPE($2)$4 ! Replaced by Perl Script(3)\n";
         next;
         }

      if ( m/^!.*Hint.*/i ||
           m/.*DIFFSIZES.*/i ) 
         {
         next;
         }

#      if ( m/^(.*)CALL(.*)_C_(DIFF|BACK)(.*)/i ) 
#         {
#         print $out "$1CALL$2$4 ! Replaced by Perl Script\n";
#         next;
#         }


      if ( m/^(.*)_C_(DIFF|BACK)(.*)/i ) 
         {
         #print "$1$3";
         print $out "$1$3 ! Replaced by Perl Script(4)\n";
         next;
         }

      if ( m/^.*POINTER.*/i ) 
         {
         next;
         }

      print $out $_;

   }

   close $out;

}

close $diffvars;

#-----------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------

foreach my $file ( <*.f90> , <*.f90> ) {

   if ( $file =~ /.*mpi_(back|diff).*/ ) {
      unlink $file
   }

   if ( $file !~ /.*_(back|diff).*/ ) {
      unlink $file
   }

   if ( $file =~ /.*_c_(back|diff).*/ ) {
      unlink $file
   }

}

#-----------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------

foreach my $file ( <*_diff.f90> , <*_back.f90> ) {

   open(my $in, "<$file");

   my @copy = <$in>;

   close $in;

   for ( @copy ) {

      if( m/.*ISIZE.*OF.*/i ) {

         print "Missing in $file : $_";

      }

   }

}
