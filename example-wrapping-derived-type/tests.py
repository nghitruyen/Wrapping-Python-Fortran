from mytype import mytype_mod

#Initialisation mytype
mine = mytype_mod.mytype(2) 

print("Valeurs initiaux mytype:")
print("a =",mine.a)
print("strickler_type =",mine.strickler_type) 
print("m =",mine.m)

print("Valeur mytype2:")
#Initialisation mytype2
mine2 = mytype_mod.mytype2(3,2) 

#Set valeurs
mytype_mod.set_m(mine2.typ[0],m=[23,5,98])

print(mine2.typ[0],mine2.typ[1])


