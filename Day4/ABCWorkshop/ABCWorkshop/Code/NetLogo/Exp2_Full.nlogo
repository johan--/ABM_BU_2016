;;;; Definition of Model Parameters ;;;;

globals [ empirical_data food_per_patch half_per_patch speed_in_patches reference_T Boltz T Arrhenius adult_size juvenile_size cocoon_size ]

patches-own 
[ change_in_food               ;;;; change in food density per timestep
  current_food                 ;;;; food density (g)
  func_response                ;;;; relative functional response
  ]

turtles-own 
[ ;;;; State Variables ;;;;
  energy_reserve              ;;;; energy reserve (kJ) : amount of energy stored as tissue (7 kJ/g)
  mass                        ;;;; individual mass (g)
  
  ;;;; Energy Management ;;;;
  energy_left                 ;;;; intermediate daily energy reserve (kJ) which is the difference between assimilated and expended energy for an individual per day. Any remaining energy enters the energy reserves.    
  energy_reserve_max          ;;;; maximum energy reserve (kJ) : a maximum threshold for the storage of energy per individual (dependant on mass)
  ingestion_rate              ;;;; ingestion rate (g/day) : amount of food ingested by an individual per day
  
  ;;;; Survival/Maintenance;;;;
  BMR                         ;;;; energy cost of maintenance (kJ) : this must be fulfilled per individual per day for survival

  ;;;; Reproduction Parameters ;;;;
  hatchlings_c                ;;;; cumulative number of hatchlings produced per adult
  hatchlings_n                ;;;; hatchlings since last measurement
  R                           ;;;; Max_R values accumulate until enough energy is available to produce one cocoon (mass_cocoon * (energy_tissue + energy_synthesis)
  ]

breed [ adults adult ]
breed [ juveniles juvenile ]
breed [ cocoons cocoon ]

to setup-interface
  set initial_number_juveniles 10
  set initial_number_adults 0
  set total_food 40
  set scape_size 0.0144
  set temperature 25
end

to setup-parameters
  set energy_food 10.6            ; kJ/g 
  set B_0 967                       ; kJ/day
  set activation_energy 0.25        ; eV
  set energy_tissue 7                ; kJ/g
  set energy_synthesis 3.6          ; kJ/g
  
  set half_saturation_coef 3.5      ; g
  set max_ingestion_rate 0.7        ; g/day 
  set mass_birth 0.011              ; g
  set mass_cocoon 0.015             ; g
  
  set mass_sexual_maturity 0.25     ; g
  set mass_maximum  0.5             ; g
  set growth_constant  0.177        ; g
  set max_reproduction_rate 0.182   ; kJ/g/day
  set speed 0.004                   ; m/day
end

to setup
  clear-all
  
  set empirical_data [ 0 30 31 39 36 41 36 48 40 44 42 38 36 42 41 34 40 36 40 ]
  
  set-default-shape adults "worm"
  set-default-shape juveniles "worm"
  set-default-shape cocoons "dot"
  
  set reference_T 298.15            ; Kelvins
  set Boltz (8.62 * (10 ^ -5))      ; eV K-1
  
  set food_per_patch total_food / count patches
  set half_per_patch ((half_saturation_coef * scape_size) / 0.01) / count patches
  set speed_in_patches speed / (sqrt(scape_size) / sqrt(count patches))
  
  set adult_size 0.9
  set juvenile_size 0.6
  set cocoon_size 0.3
  setup-patches
  setup-turtles
  
  set T 273.15 + temperature
  set Arrhenius (e ^ ((- activation_energy / Boltz ) * ((1 /  T ) - (1 / reference_T))))
  
  reset-ticks
end

to setup-patches
  ;;; each patch is asked to set its colour as green and its food density to that of 'food_density_patch' the value of which is determined by the user on the interface slider
  ask patches 
  [ set current_food food_per_patch

    ifelse current_food > 0 
    [ set pcolor scale-color green current_food (food_per_patch * 2) (0 - food_per_patch * 0.5) ]
    [ set pcolor brown ]
  
    calc-func-response ]
end

to setup-turtles
  ;;; the initial population density of each life cycle stage is determined by the user on the interface
  ;;; colour, size and state variables are set depending on the life cycle stage
  ;;; all individuals are set at a random position within the landscape
  
  create-adults initial_number_adults
  [ set color red
    set size adult_size
    set mass mass_sexual_maturity
    setxy random-xcor random-ycor ]
  
  create-juveniles initial_number_juveniles
  [ set color pink
    set size juvenile_size
    set mass mass_birth
    setxy random-xcor random-ycor ]

   ask turtles with [ breed != cocoons ]
   [ set energy_reserve (mass / 2) * energy_tissue ]
end
  
to go ;;; when the go button on the interface is pressed, the following schedule of processes occurs in one timestep.
  if not any? turtles with [ breed != cocoons ] or ticks > 210 [ stop ]

  if ticks = 25
  [ ask adults
    [ set hatchlings_c hatchlings_c + hatchlings_n
      set hatchlings_n 0 ]
    ask cocoons [ die ] ]
  
  if ticks > 25 and remainder ticks 10 = 5
  [ ask adults
    [ set hatchlings_c hatchlings_c + hatchlings_n
      set hatchlings_n 0 ]
    ask cocoons [ die ] ]
  
  ask turtles with [ breed != cocoons ]
  [ set energy_reserve_max (mass / 2) * energy_tissue
    calc-ingestion-rate ]
    
  ask patches [ correct-ingestion-rate ]
  
  ask adults
  [ calc-assimilation
    calc-maintenance
    calc-reproduction
    calc-growth
    update-reserves
    move ]
  
  ask juveniles
  [ calc-assimilation
    calc-maintenance
    calc-growth
    transform-juvenile
    update-reserves
    move ]
  
  ask patches
  [ update-patch ]
  
  tick
end

to go-25 
  repeat 25 [ go ]
end

to go-10 
  repeat 10 [ go ]
end

;;;;;;;;;;;;;;;;;; Ingestion Rate ;;;;;;;;;;;;;;;;;;;;
;;; juveniles and adults calculate their ingestion rate (the amount of food ingested from the environment) which depends on the food density of the patch in which they are present and the mass dependent maximum ingestion rate.
to calc-ingestion-rate
   set ingestion_rate (max_ingestion_rate * Arrhenius) * func_response * (mass ^ (2 / 3))
end

;;;;;;;;; Assimilation of Energy from Ingested Food ;;;;;;;;;
;;; assimilation is the amount of energy available from the given amount of ingested food, determined by the energy content of food and assimilation efficiency, irrespective of mass or temperature 
;;; Energy_assimilated provides an intermediate energy reserve, where the energy assimilated in one time step may be used by metabolic processes before being fixed within the energy reserves

to calc-assimilation
  set energy_left (ingestion_rate * energy_food)
end

;;;;;;;;;;;;;;;;;;;;; Somatic Maintenance ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; maintenance costs are calculated as BMR, which is essential for individual survival and is mass and temperature dependant
to calc-maintenance
  set BMR B_0 * (mass ^ (3 / 4)) * e ^ (- activation_energy / (Boltz * T))     ; { Equation 1 }
  
  ;;; if energy assimilated in one time step is available to pay for maintenance, the costs are deducted from this intermediate reserve. If energy is not available the energy costs are depleted from the reserves in order to maximise survival
  ifelse energy_left > BMR
  [ set energy_left energy_left - BMR ]
  [ set energy_reserve energy_reserve - (BMR - energy_left)
    set energy_left 0 ]
  
  ;;; if the conditions are such that scaled energy reserves fall below the scaled mass of an individual that individual is said to be starving and is asked to undergo processes to offset mortality through a starvation strategy
  if energy_reserve < (energy_reserve_max * 0.5) and breed != cocoons
  [ onset-starvation-strategy ]
end
       
;;;;;;;;;;;;;;;;;;;; Starvation Strategy ;;;;;;;;;;;;;;;;;;;;;;;;;
to onset-starvation-strategy
  ;;; the starvation strategy asks individuals to lose weight in order to catabolise energy for covering maintenance costs. This is taken from the mass of the individual and the energy content of the lost mass becomes available within the energy reserves to cover maintenance costs
  ifelse (energy_tissue + energy_synthesis > 0)
  [ set mass mass - (BMR / (energy_tissue + energy_synthesis)) ]
  [ die ]
  set energy_reserve energy_reserve + BMR
  
  ;;; if the mass of an individual adult falls below the mass at puberty due to starvation these individuals rejuvenate to their pre-sexual mature state and are classified as juveniles, unable to reproduce
  if mass < mass_sexual_maturity 
  [ set breed juveniles
    set color pink
    set size juvenile_size ]
  
  ;;; if the mass of any individual falls below the mass at birth they die due to their inability to sustain such weight loss 
  if mass < mass_birth 
  [ die ]  
end

;;;;;;;;;;;;;;;;;;; Reproduction ;;;;;;;;;;;;;;;;;;;;
;;; Energy available after maintenance in adults goes to ova development. Energy accumulates until enough is available to produce one cocoon containing one fully developed ova

to calc-reproduction
  let max_R (max_reproduction_rate * Arrhenius) * mass
  
  let energy_for_R min list energy_left max_R
  set energy_left energy_left - energy_for_R

  if (energy_for_R < max_R and energy_reserve > 0.5 * energy_reserve_max)
  [ set energy_reserve energy_reserve - (max_R - energy_for_R)
    set energy_for_R max_R ] 
  
  set R R + energy_for_R
  
  if R >= (mass_cocoon * (energy_tissue + energy_synthesis))
  [ reproduce ]
end

to reproduce
  hatch-cocoons 1
  [ set color white
    set size cocoon_size ]
  
  set hatchlings_c hatchlings_c + 1
  set hatchlings_n hatchlings_n + 1
  set R (R - (mass_cocoon * (energy_tissue + energy_synthesis))) 
end

;;;;;;;;;;;;;;;;;;;;; GROWTH ;;;;;;;;;;;;;;;;;;;;;;;;
;;; after maintenance in juveniles and maintenance and reproduction in adults, available energy is expended on growth
;;; a maximum growth rate, following the von Bertalanffy growth equation 
to calc-growth
  let max_G (growth_constant * Arrhenius) * (mass_maximum ^ (1 / 3) * mass ^ (2 / 3) - mass)
  let energy_for_G min list (max_G * (energy_tissue + energy_synthesis)) energy_left
  
  let to_grow 0
  
  if (max_G > 0 and (energy_tissue + energy_synthesis) > 0)
  [ set to_grow (energy_for_G / (max_G * (energy_tissue + energy_synthesis))) * max_G ]  
  
  if (mass + to_grow) < mass_maximum
  [ set mass mass + to_grow
    set energy_left energy_left - energy_for_G ]
end

to update-reserves
  ;;; energy remaining after metabolic costs have been covered are fixed within the energy reserves, with the energy costs of synthesis taken into account as energy is stored as flesh
  
  if energy_left > 0
  [ ifelse (energy_tissue + energy_synthesis > 0)
    [ set energy_reserve energy_reserve + (energy_left * (energy_tissue / (energy_tissue + energy_synthesis))) ]
    [ die ] ]
  
  ;;; a maximum threshold of 50% mass * energy content of flesh is set for energy reserves   
  if energy_reserve > energy_reserve_max and breed != cocoons 
  [ set energy_reserve energy_reserve_max ]
end

;;;;;;;;;;;;;;; Movement of Individuals ;;;;;;;;;;;;;;;;;;;;;;;;;
;;; juveniles and adults are asked to move in a random direction across the landscape
to move   
    rt random 90
    lt random 90
    fd speed_in_patches
end

;;;;;;;;;;;;;;;;; Life Stage Transformations ;;;;;;;;;;;;;;
;;; cocoons transform to the juvenile stage when their age is equivalent to the temperature dependent incubation period and embryonic development is 100%

;;; juveniles transform to the adult life stage when they have grown to a mass equivalent to that at puberty (mass_sexual_maturity)
to transform-juvenile
  if mass >= mass_sexual_maturity 
  [ set breed adults
    set color red
    set size adult_size ]
end

;;;;;;;;;;;;;;;;; Updating a Patch ;;;;;;;;;;;;;;
;; at the end of every time step, the food, functional response and color of patches is updated
to update-patch
  calc-change-in-food
  set current_food max list (current_food - change_in_food) 0
  
  if ticks = 25 [ set current_food (20 / count patches) ]
  if ticks > 25 and remainder (ticks - 25) 20 = 0
  [ set current_food (current_food + (20 / count patches)) ]

  calc-func-response
  
  ifelse current_food > 0 
  [ set pcolor scale-color green current_food (food_per_patch * 2) (0 - food_per_patch * 0.5) ]
  [ set pcolor brown ]
end 

;;;;;;;;;;;;;;;; Change in Food Density ;;;;;;;;;;;;;;;;
;;; the change in food density per time step is equivalent to the amount individuals ingest from one patch
to calc-change-in-food
  ifelse current_food <= 0 
  [ set change_in_food 0 ]
  [ set change_in_food sum [ ingestion_rate ] of turtles-here with [ breed != cocoons ] ]
end

;;;;;;;;;;;;;;; Relative Functional Response ;;;;;;;;;;;;;;;;;
;;; the functional response is calculated for each patch so that the food ingested by individuals feeding there depends on food density
to calc-func-response
  ifelse (current_food + half_per_patch) > 0
  [ set func_response current_food / (current_food + half_per_patch) ]
  [ set func_response current_food / (current_food + 1E-10) ]
end

;;;;;;;;;;;;;;; Correcting Ingestion Rate ;;;;;;;;;;;;;;;;;
;;; if there is not enough food left for all turtles, the patch distributes its remaining food equally
to correct-ingestion-rate
  ;;; if no food is available within a patch the ingestion rate is set to 0, to avoid negative values occuring from the regression equation above
  if current_food = 0
  [ ask turtles-here with [ breed != cocoons ] [ set ingestion_rate 0 ] ]
  
  if sum [ ingestion_rate ] of turtles-here with [ breed != cocoons ] > current_food
  [ ask turtles-here with [ breed != cocoons ] [ set ingestion_rate current_food / count turtles-here with [ breed != cocoons ] ]
    set current_food 0 ]
end

;;;;;;;;;;;;;;;;; Plotting Functions ;;;;;;;;;;;;;;
;; these functions plot the lines in the display, to give an overview of the model's performance while it is running
to-report mean-mass
  ifelse any? turtles with [ breed != cocoons ]
  [ report mean [ mass ] of turtles with [ breed != cocoons ] ]
  [ report 0 ]
end

to-report sum-hatchlings
  ifelse any? turtles with [ breed != cocoons ]
  [ report sum [ hatchlings_n ] of turtles with [ breed != cocoons ] ]
  [ report 0 ]
end
@#$#@#$#@
GRAPHICS-WINDOW
256
22
566
353
-1
-1
150.0
1
10
1
1
1
0
1
1
1
0
1
0
1
1
1
1
ticks
30.0

BUTTON
594
351
658
384
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
903
352
966
385
NIL
Go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
48
26
229
59
initial_number_juveniles
initial_number_juveniles
0
10
10
1
1
NIL
HORIZONTAL

SLIDER
47
69
228
102
initial_number_adults
initial_number_adults
0
50
0
10
1
NIL
HORIZONTAL

SLIDER
44
181
228
214
temperature
temperature
0
30
25
1
1
NIL
HORIZONTAL

INPUTBOX
46
111
129
171
total_food
40
1
0
Number

PLOT
594
25
911
173
Food Available
time (days)
total food (g)
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"model" 1.0 0 -9276814 true "" "plot sum [ current_food ] of patches"

PLOT
593
187
1037
337
Number of Cocoons
time (ten days)
cocoons
0.0
10.0
0.0
50.0
true
true
"" "if ( ticks >= 25 and remainder ticks 10 = 5 )\n[ set-current-plot-pen \"model  \"\n  plot sum-hatchlings\n  set-current-plot-pen \"data\"\n  plot item ((ticks - 25) / 10) empirical_data ]"
PENS
"data" 1.0 0 -16777216 true "" ""
"model  " 1.0 0 -7500403 true "" ""

INPUTBOX
888
398
976
458
energy_food
10.6
1
0
Number

INPUTBOX
592
398
664
458
B_0
967
1
0
Number

INPUTBOX
673
398
787
458
activation_energy
0.25
1
0
Number

INPUTBOX
795
398
880
458
energy_tissue
7
1
0
Number

INPUTBOX
985
397
1090
457
energy_synthesis
3.6
1
0
Number

INPUTBOX
919
471
1002
531
mass_cocoon
0.015
1
0
Number

INPUTBOX
718
470
828
530
max_ingestion_rate
0.7
1
0
Number

INPUTBOX
591
469
709
529
half_saturation_coef
3.5
1
0
Number

INPUTBOX
837
470
911
530
mass_birth
0.011
1
0
Number

INPUTBOX
589
541
720
601
mass_sexual_maturity
0.25
1
0
Number

INPUTBOX
730
542
834
602
mass_maximum
0.5
1
0
Number

INPUTBOX
846
543
948
603
growth_constant
0.177
1
0
Number

INPUTBOX
958
543
1087
603
max_reproduction_rate
0.182
1
0
Number

BUTTON
667
351
752
384
Go Once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
826
352
895
385
Go 10
go-10
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
973
352
1115
385
Set Basic Interface
setup-interface
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1123
352
1277
385
Set Basic Parameters
setup-parameters
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
759
351
819
384
Go 25
go-25
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1097
543
1182
603
speed
0.0040
1
0
Number

INPUTBOX
139
111
228
171
scape_size
0.0144
1
0
Number

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

worm
true
0
Polygon -7500403 true true 165 210 165 225 135 255 105 270 90 270 75 255 75 240 90 210 120 195 135 165 165 135 165 105 150 75 150 60 135 60 120 45 120 30 135 15 150 15 180 30 180 45 195 45 210 60 225 105 225 135 210 150 210 165 195 195 180 210

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 1.0 0.0
0.0 1 1.0 0.0
0.2 0 1.0 0.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
