# Wrapping-Python-Fortran

Ce projet est une partie des missions de mon stage chez CNRS à partir de Juin à Septembre 2020 avec le sujet "Wrapping du code de DassFlow2D".
Pour information, le DassFlow2D est un logiciel écrit principalement en Fortran, dédié à la simulation hydraulique fluviale et notamment, conçu pour l'assimilation variationnelle des données.
Le but du stage est de faire une conversion du code de Python-Fortran (on appelle le wrapping) du code de DassFlow2D pour que l'on puisse appeler et utiliser des routines et des variables en Fortran afin d’effectuer une simulation simple en Python.

Ce projet "Wrapping Python-Fortran" contient des exemples simples pour tester la fonctionnalité du code de wrapping avec les générateurs f2py et notamment f90wrap.
1. Le répertoire "f90wrap-test" contient des cas de test avec le niveau de complexité augmente de /ex_1 à /ex_3.
2. Le répertoire "wrapping-test-with-makefile" est un test pour faire tourner automatiquement le code de wrapping avec un Makefile.
3. La fonctionnalité du cas de type dérivé avec le Makefile est dans /example-wrapping-derived-type.
4. Le "dassflow2d-wrapped" contient un extrait des fichiers importants pour le wrapping de DassFlow2D et le fichier Python pour faire une simulation.
