globals
[collisions Mean_tr_time Sum_tr_time amount_gone_cells Mean_sch_time Sum_sch_time amount_cont_cells am_gone_red_cells dying_cells
  x_free y_free    ; coordinates of neighbour non-occupied patch
 ;thresh_level     ; activation level's treshold value, after which division can occur
  tn               ; Normalized contact time
 ;mean_aff         ; Mean affinity value (exponentially distributed) for T cell population
  ]

turtles-own
[ angle contacts cnt_1 cnt_3 cnt_4 cnt_5 cnd
  contact_flag    ; Was or not contact with Dendritic cell for this T cell [0, 1]
  tr_time         ; Transit time through computational domain
  crt_time        ; Time of lymphocyte occurence in the system
  sch_time        ; T cell - DC search time
  affinity        ; Cognate or non-cognate lymphocyte: float from 0 to 1
  activity        ; Activate or non-activate T cell
  effectiveness   ; Can T cell divide or already effector [0, 1]
  doubling_time   ; Division time
  nn              ; Number of turtles-neighbours
  nbp             ; Number of neighbour non-occupied patches
  nh              ; Number of turtles on the same patch
  num_gen         ; Generation number for dividing T cells
  max_div         ; Divisions number limit
  state           ; 0 -- free state, outside gradient
                  ; 1 -- DC contact
                  ; 3 -- desensitization state
                  ; 4 -- gradient-sensitive state
                  ; 5 -- waiting for division

  activation_time ; T cell priming time
  contact_time    ; Duration of short contacts for cognate T cells
  Prime_prob      ; Probability for T cell to be primed as a function of affinity
  Prime_flag      ; Will a naive cell be primarily activated during this contact? (0 - no, 1 - yes)
  a_level         ; T cell activation level
  zero_level      ; Initial T cell activation level to the moment of DC contact start
  TB              ; Lognormally distributed random times of short contacts independent from affinity
  ]

patches-own
[ sinusity
  DC_center       ; = 1 for DC center
  stop-zone       ; The zone in which the centers of other DCs cannot be located, so that they are not located too close to each other.

 ]


to setup

  clear-all

  setup-patches

  setup-turtles

  draw-walls

  reset-ticks

end
;
;
;
;
to setup-patches                 ; DC construction
  ask patches [set pcolor blue]
  let i 0
  let j 0
  repeat 200
  [
    ask patch (-44 + random 88 ) (-44 + random 88 )
    [ if( stop-zone != 1 )
      [ set i i + 1
        let px pxcor
        let py pycor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ifelse show-gradient?
        [
          ask patch px py
          [ ask patches in-radius 5
             [ set pcolor 4 ]
            ask patches in-radius 20
             [ set stop-zone 1 ]
          ]
        ]
        [
          ask patch px py
          [ ask patches in-radius 5
            [ set pcolor blue ]
            ask patches in-radius 20
             [ set stop-zone 1 ]
          ]
        ]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ask patch px py
        [ set pcolor black
          set DC_center 1  ]

        ask patch (px - 1) py
        [ set pcolor black ]

        ask patch (px + 1) py
        [ set pcolor black ]

        ask patch px (py - 1)
        [ set pcolor black ]

        ask patch px (py + 1)
        [ set pcolor black ]

      ]]
    if (i = 8) [stop]
  ]
end
;
;
;
;
;
;
;
to draw-walls
  ask patches with [abs pxcor = 50] [ set pcolor violet ] ; designation of side boundaries with periodic conditions
  ask patches with [pycor = 50] [ set pcolor green ]  ; reflective wall of the capsule was drawn
  ask patches with [pycor = -50] [ set pcolor green ] ; reflective wall of the capsule was drawn

  if what-a-boundary? = "entire adsorbing border"
  [  ask patches with [abs pycor = 50]
      [ set pcolor blue
        set sinusity 1
      ] ; all bottom boundary is free for T cell's exit
  ]

  if what-a-boundary? = "120 patches"
  [  ask patches with [abs pycor = 50] with [pxcor > -40]  with [pxcor < -10]
      [ set pcolor blue
        set sinusity 1
      ]

     ask patches with [abs pycor = 50] with [pxcor > 10]  with [pxcor < 40]
     [ set pcolor blue
       set sinusity 1
     ]
  ]

  if what-a-boundary? = "28 patches"
  [  ask patches with [abs pycor = 50] with [pxcor > -27]  with [pxcor < -20]
      [ set pcolor blue
        set sinusity 1
      ]

     ask patches with [abs pycor = 50] with [pxcor > 20]  with [pxcor < 27]
     [ set pcolor blue
       set sinusity 1
     ]
  ]

  if what-a-boundary? = "20 patches"
  [  ask patches with [abs pycor = 50] with [pxcor > -30]  with [pxcor < -25]
      [ set pcolor blue
        set sinusity 1
      ]

     ask patches with [abs pycor = 50] with [pxcor > 25]  with [pxcor < 30]
     [ set pcolor blue
       set sinusity 1
     ]
  ]

    if what-a-boundary? = "60 patches"
  [  ask patches with [abs pycor = 50] with [pxcor > -35]  with [pxcor < -20]
      [ set pcolor blue
        set sinusity 1
      ] ; Medullary sinus was drawn

     ask patches with [abs pycor = 50] with [pxcor > 20]  with [pxcor < 35]
     [ set pcolor blue
       set sinusity 1
     ]
  ]

  if what-a-boundary? = "52 patches"
  [  ask patches with [abs pycor = 50] with [pxcor > -33]  with [pxcor < -20]
      [ set pcolor blue
        set sinusity 1
      ] ; Medullary sinus was drawn

     ask patches with [abs pycor = 50] with [pxcor > 20]  with [pxcor < 33]
     [ set pcolor blue
       set sinusity 1
     ]
  ]

      if what-a-boundary? = "64 patches"
  [  ask patches with [abs pycor = 50] with [pxcor > -36]  with [pxcor < -20]
      [ set pcolor blue
        set sinusity 1
      ] ; Medullary sinus was drawn

     ask patches with [abs pycor = 50] with [pxcor > 20]  with [pxcor < 36]
     [ set pcolor blue
       set sinusity 1
     ]
  ]

      if what-a-boundary? = "40 patches"
  [  ask patches with [abs pycor = 50] with [pxcor > -35]  with [pxcor < -25]
      [ set pcolor blue
        set sinusity 1
      ] ; Medullary sinus was drawn

     ask patches with [abs pycor = 50] with [pxcor > 25]  with [pxcor < 35]
     [ set pcolor blue
       set sinusity 1
     ]
  ]

      if what-a-boundary? = "32 patches"
  [  ask patches with [abs pycor = 50] with [pxcor > -34]  with [pxcor < -26]
      [ set pcolor blue
        set sinusity 1
      ] ; Medullary sinus was drawn

     ask patches with [abs pycor = 50] with [pxcor > 26]  with [pxcor < 34]
     [ set pcolor blue
       set sinusity 1
     ]
  ]

  if what-a-boundary? = "8 patches"
  [  ask patches with [abs pycor = 50] with [pxcor > -2]  with [pxcor < 2]
      [ set pcolor blue
        set sinusity 1
      ] ; Medullary sinus was drawn
  ]

    if what-a-boundary? = "no passage"
  [  ask patches with [abs pycor = 50]
      [ set pcolor green
        set sinusity 0
      ]                             ; If there is no exit for T cells
  ]
end
;
;
;
;

to setup-turtles
  set-default-shape turtles "circle"
  ;set mean_aff 0.5
  crt number + 7
  [
   set color grey
   set activity 0      ;
   set effectiveness 0 ;
   set state 0
   set contact_flag 0
   set Prime_flag 0
   set crt_time 0
   set size 1
   set num_gen 0
   setxy random-xcor random-ycor

   set doubling_time round random-normal 960 150
   set activation_time round random-normal 2880 250
   set max_div 2 + random 17

   let A random-exponential mean_aff
   ifelse A > 1
      [ set affinity 1 ] ; instead of random-float 1
      [ set affinity A ]

  ]

  ask turtles [
   if
   pcolor = black [die]
   ]
end
;
;
;
;
to
  go
  ; Stop after 120000 ticks
  if
  ticks >= sim_time [ stop ]

  move-nonactivate
  move-turtles
  move-effector

  gradient-move

  free-move-turtles

  tick

  count-collisions

  chemokine_overdose
  chemokine_reset
  saturation-count

  contact-change
  contact-count

  wait-change
;  prime
;  free-lymphocyte
;  wait-count

  change-state
  duration-count

  wake

  exceed-threshold
  wait-replication
  divide

  life-count
  kill-effectors

  recolor

  calculate_cognate_cells

  calculate-level

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Part about T cell movement
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to
  move-turtles
  ;
  ask turtles with [pcolor != 4] with [activity = 1]
  [

    rt random 80
    lt random 80

    if not any? (turtles-on patch-ahead 1) with [self != myself] and [pcolor] of patch-ahead 1 != black and [pcolor] of patch-ahead 1 != green and state != 1
      [ fd 1 ]
  ]
   turtles-out
end

to
  move-nonactivate ; modified

  ask turtles with [ activity = 0 ]
  [

    rt random 80
    lt random 80

    if not any? (turtles-on patch-ahead 1) with [self != myself] and [pcolor] of patch-ahead 1 != black and [pcolor] of patch-ahead 1 != green and state != 1
      [ fd 1 ]
  ]

   turtles-out
end


to
  move-effector
  ;
  ask turtles with [ effectiveness = 1]
  [

    rt random 80
    lt random 80

    if not any? (turtles-on patch-ahead 1) with [self != myself] and [pcolor] of patch-ahead 1 != black and [pcolor] of patch-ahead 1 != green and state != 1
      [ fd 1 ]
  ]

   turtles-out
end


to gradient-move               ; modified
                               ; non-random movement in the chemokine clouds

  ask turtles with [pcolor = 4 ] with [activity = 1] with [ (state = 4 or state = 5 )] with [effectiveness = 0]
  [ ifelse random Prob_chem = 0

      [ let target-patch min-one-of (patches with [DC_center = 1]) [distance myself]
        if target-patch != nobody
          [ set heading towards target-patch ]
        if not any? (turtles-on patch-ahead 1) with [self != myself] and [pcolor] of patch-ahead 1 != black
          [fd 1]
      ]

      [
        rt random 80
        lt random 80
        if not any? (turtles-on patch-ahead 1) with [self != myself] and [pcolor] of patch-ahead 1 != black
         [fd 1]
      ]
   ]
  ;
end


;
;
;
;
to
  free-move-turtles              ; свободное движение в состоянии десенситизации

  ask turtles with [pcolor = 4] [

    rt random 80
    lt random 80

    if not any? (turtles-on patch-ahead 1) with [self != myself] and [pcolor] of patch-ahead 1 != black and [pcolor] of patch-ahead 1 != green and (state = 3)
      [ fd 1 ]
    ]
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to measure_tr_time ; Transit time measurements
  ask turtle 1
  [
    if ( [pycor] of patch-ahead 1 = -50)
    [
      set tr_time ticks
      set Mean_tr_time Mean_tr_time + tr_time
    ]
  ]
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;  Part about T cell exit through MS and new cells occurence
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to
  turtles-out                                    ; T-cells go to the bloodstream through the lower boundary, and immediately
                                                 ; there are new T-cells from the random place of the Computational domain only if cell count < 2000 cells

  ask turtles
  [
    if ( [sinusity] of patch-ahead 1 = 1)
    [ if count turtles < pop_limit
      [
        ask one-of patches with [pcolor = blue] with [count turtles-here = 0]           ; patch -50 random-ycor
        [
          sprout 1
            [ set color grey

              set activity 0
              set num_gen 0
              set Prime_flag 0
              set effectiveness 0
              set size 1
              rt random-float 360
              set crt_time ticks
              set doubling_time round random-normal 960 150
              set activation_time round random-normal 2880 250
              set max_div 2 + random 17

              let A random-exponential mean_aff
              ifelse A > 1
              [ set affinity 1 ]
              [ set affinity A ]

            ]
          ]
      ]
      set tr_time ticks - crt_time
      set amount_gone_cells amount_gone_cells + 1

      set Sum_tr_time Sum_tr_time + tr_time

      if what-a-output? = "Transit time"
      [output-print tr_time]

      if what-a-output? = "Cognate clones output"
        [ output-type ticks / 2
          output-type " "
          output-print am_gone_red_cells
        ]

      if what-a-output? = "Affinity distribution"
      [
        if activity = 1
          [output-print affinity]
      ]

      if what-a-output? = "Affinity distribution effector"
      [
        if effectiveness = 1
          [output-print affinity]
      ]
      die
      ]
    ]
end


to count-collisions ; Count contacts and unique contacts of T cells with DCs
  ask turtles
  [

   if [pcolor] of patch-ahead 1 = black and contact_flag = 0

     [ set collisions collisions + 1
       set contacts contacts + 1
       set sch_time ticks - crt_time
       set Sum_sch_time Sum_sch_time + sch_time
       set amount_cont_cells amount_cont_cells + 1
       if what-a-output? = "Unique contacts"
       [ output-type ticks / 2
         output-type " "
         output-print collisions
       ]
     ]
   ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Part for T cell's state changes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to chemokine_overdose                                            ; The start of the countdown to the loss of the chemokine sensitivity due to too long stay in the chemokine cloud
  ask turtles with [activity = 1] with [effectiveness != 1]
  [ if pcolor = 4 and state = 0
    [
      set state 4
      ;set color orange + 3
    ]
  ]
end

to chemokine_reset                                          ; Reset the state of desensitizations waiting and the corresponding counter to zero
                                                            ; after the lymphocyte has left the gradient cloud with a counter filled less than half
  ask turtles with [activity = 1] with [effectiveness != 1]
  [ if pcolor = blue and state = 4 and cnt_4 <= 10          ; It is hardly necessary to continue counting the waiting for desensitization, if the lymphocyte is no longer in a gradient.
    [
      set state 0
      set cnt_4 0
    ]
  ]
end

to
  saturation-count
  ask turtles with [activity = 1] with [effectiveness != 1]
   [ if state = 4
     [
       set cnt_4 cnt_4 + 1
     ]
   ]
end

to
  contact-change                                         ; Come into DC contact state

  ask turtles with [effectiveness = 1]
  [ if [pcolor] of patch-ahead 1 = black and state != 1
    [
      set state 1
      set contact_flag 1
    ]
  ]


  ask turtles with [effectiveness != 1] with [activity = 0]
  [
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      ; Hill function dependency from T cell's affinity for priming probability
       let n 5
       set Prime_prob (affinity ^ n) / (affinity ^ n + K_prime ^ n)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    if [pcolor] of patch-ahead 1 = black and state != 1 and state != 5
      [
       set state 1
       set contact_flag 1
       set zero_level a_level
       ifelse Prime_prob > random-float 1
        [ set Prime_flag 1 ]
        [ set Prime_flag 0 ]
      ]
  ]

  ask turtles with [effectiveness != 1] with [activity = 1]
  [ if [pcolor] of patch-ahead 1 = black and state != 1 and state != 5
      [set state 1
       set contact_flag 1
       set zero_level a_level
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      ; Hill function dependency from T cell's affinity for contact time and priming probability
         let n 5
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      ; 20 min and dispersion for log-normal distribution construction
         let Mu 100
         let sigma 50
         let beta ln(1 + sigma * sigma / (Mu * Mu))
         let M ln Mu - beta / 2
         let S sqrt beta
         set TB round exp(random-normal M S)
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       set tn (affinity ^ n) / (affinity ^ n + K_time ^ n)
       set contact_time ceiling (40 * tn)
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ]
  ]
end

to
  contact-count
  ask turtles
   [ if state = 1
     [
       set cnt_1 cnt_1 + 1
     ]
   ]
end

to
  wait-change ; What happens after the contact with the DC was finished? Transition to baseline state, waiting for division or desensitization.

  ask turtles with [effectiveness = 1] ; For effectors
  [
    if cnt_1 >= 6
      [
        set state 0
        set cnt_1 0
      ]
  ]

  ask turtles with [activity = 0] with [effectiveness != 1] with [Prime_flag = 0]; For non-activated naive T cells, which will not be primed during this contact
  [
    if cnt_1 >= 6
      [
        set state 0
        set cnt_1 0
      ]
  ]

  ask turtles with [activity = 0] with [effectiveness != 1] with [Prime_flag = 1]; For non-activated naive T cells, which will be primed during this contact
  [
    if cnt_1 >= activation_time        ; (Duration of priming is 24 hours.)
     [
          set state 5
          set cnt_1 0
          set activity 1
          set Prime_flag 0
          set color red
      ]
  ]

  ask turtles with [activity = 1] with [effectiveness != 1] ; For already activated T cells.
  [
    if ( cnt_1 >= contact_time * CTLA4_time )                                   ; "Short" contacts of already activated T cells last a random time and determine the level of activation.
    [
      ifelse ( a_level > thresh_level )
        [
          set state 5
          set cnt_1 0
        ]
        [
          set state 3
          set cnt_1 0
        ]
    ]
  ]
end

to
  change-state ; Transition to the desesitization state

  ask turtles with [activity = 1] with [effectiveness != 1] with [state = 4]
    [
     if cnt_4 = 20
       [
        set state 3
        set cnt_4 0
       ]
    ]
end

to
  duration-count  ; Duration of desensitization
  ask turtles with [activity = 1] with [effectiveness != 1]
   [ if state = 3 and pcolor = blue
     [
       set cnt_3 cnt_3 + 1
     ]
   ]
end

to
  wake  ; Transition from a desensitization chemokine state to a base chemokine-sensitive state
  ask turtles with [activity = 1] with [effectiveness != 1]
    [
     if cnt_3 = 20
       [
        set state 0
        set cnt_3 0

       ]
    ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Code part on the cognate cells division
to
  calculate-neighbors
   ask turtles
     [
      set nn count turtles-on neighbors

      set nh count turtles-on patch-here


      ;set nbp count neighbors with [pcolor = black]
     ]

end

to exceed-threshold    ; accelerated transition to the state of division waiting for all cells, whose a_level > threshhold level
  ask turtles with [activity = 1] with [effectiveness != 1] with [state != 5] with [state != 1]
    [
      if a_level >= thresh_level
        [
          set state 5
          set cnt_1 0
        ]
    ]
end

to wait-replication  ; Time to the division start
  ask turtles
   [ if state = 5
     [
       set cnt_5 cnt_5 + 1
     ]
   ]
end

to
  divide  ; After 8 hours, division occurs.
  ask turtles with [state = 5] with [num_gen <= max_div - 1 ]
    [

;;;;;;;;;;;;;;;;;;;;;;;;;;;;Search-place

       let free-patches patches in-radius Division_rad with [pcolor != black] with [count turtles-here = 0]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; replicate

     if cnt_5 >= doubling_time and any? free-patches
       [set state 0
        set cnt_5 0
        set contact_flag 1
        set activity 1
        set num_gen num_gen + 1

        hatch 1
         [ ;setxy x_free y_free
           set color red
           set affinity affinity
           set activity 1                                      ;;;;;;;;;;;;;;; New activated T cell occurs
           set effectiveness 0
           set num_gen num_gen
           set state 0
           set cnt_5 0
           set cnd 0
           set size 1
           set crt_time ticks
           set contact_flag 0
           set Prime_flag 0
           set doubling_time round random-normal 960 150                   ; New activated T cell occurs => you must immediately set the division time
           set activation_time round random-normal 2880 250
           set max_div 2 + random 17
           set a_level a_level
           move-to one-of free-patches ;;;
         ]

       ]
    ]
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Code part for life-time count

to
  life-count
  ask turtles with [effectiveness = 1]
   [
     set cnd cnd + 1
   ]
end

to
  kill-effectors
  ask turtles with [effectiveness = 1]
   [
     if cnd > life-time
       [die]
   ]
end

to recolor
  ask turtles with [activity = 1]
    [
     if num_gen >= max_div
       [
         set color orange
         set effectiveness 1
       ]

    ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; End of T cell state's changig code part
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calculate-level                                                        ; Calculation of T cell activation level

  ask turtles
   [
     ifelse state = 1
     [set a_level zero_level + alfa / ( 1 + exp (- beta_2 * cnt_1) ) ]

     [set a_level a_level - 0.000241 * lambda * a_level]                  ; Decay rate for half-life = 24 h multiplied by lambda
   ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to
  calculate_cognate_cells
    if what-a-output? = "Cognate clones number"
      [ output-type ticks / 2
        output-type " "
        output-print count turtles with [activity = 1]
      ]
end
@#$#@#$#@
GRAPHICS-WINDOW
401
84
1015
699
-1
-1
6.0
1
10
1
1
1
0
1
1
1
-50
50
-50
50
1
1
1
ticks
30.0

BUTTON
76
79
139
112
NIL
setup
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
76
132
139
165
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1074
121
1200
166
Activated cells
count turtles with [color = red]
17
1
11

PLOT
1131
313
1601
706
plot 1
Time
Number of cognate cells
0.0
20.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -2674135 true "" "plot count turtles with [color = red]"
"pen-1" 1.0 0 -955883 true "" "plot count turtles with [color = orange]"
"pen-2" 1.0 0 -16777216 true "" "plot count turtles with [activity = 1]"

INPUTBOX
178
310
286
370
number
2000.0
1
0
Number

MONITOR
1197
42
1316
87
NIL
amount_gone_cells
17
1
11

MONITOR
1437
45
1553
90
NIL
amount_cont_cells
17
1
11

SWITCH
42
557
183
590
show-gradient?
show-gradient?
0
1
-1000

CHOOSER
4
603
183
648
what-a-boundary?
what-a-boundary?
"entire adsorbing border" "120 patches" "64 patches" "60 patches" "52 patches" "40 patches" "32 patches" "28 patches" "20 patches" "8 patches" "no passage"
6

INPUTBOX
179
106
290
166
sim_time
80640.0
1
0
Number

TEXTBOX
179
87
310
105
Simulation time (in ticks):
11
0.0
1

OUTPUT
1667
105
1860
708
12

CHOOSER
1668
45
1870
90
what-a-output?
what-a-output?
"Affinity distribution" "Affinity distribution effector" "Transit time" "Unique contacts" "Cognate clones output" "Cognate clones number" 3
1

MONITOR
1442
128
1624
173
Cognate clones output:
am_gone_red_cells
17
1
11

INPUTBOX
181
409
290
469
Prob_chem
3.0
1
0
Number

TEXTBOX
176
375
326
403
Probability of chemotaxis driven motion towards DC:
11
0.0
1

TEXTBOX
178
292
299
310
Initial number of T cells:
11
0.0
1

MONITOR
1442
204
1630
249
Died cognate clones number:
dying_cells
17
1
11

INPUTBOX
25
410
130
470
Division_rad
3.0
1
0
Number

TEXTBOX
23
377
160
419
Search radius of free patch during cell division:
11
0.0
1

INPUTBOX
25
224
131
286
lambda
3.0
1
0
Number

TEXTBOX
38
193
130
235
Stimulation signal decay rate:
11
0.0
1

SLIDER
42
502
183
535
thresh_level
thresh_level
0
10
2.0
0.1
1
NIL
HORIZONTAL

INPUTBOX
23
311
129
371
life-time
7200.0
1
0
Number

TEXTBOX
35
293
128
311
Effector life-time:
11
0.0
1

MONITOR
642
16
761
73
T cells in LN
count turtles
17
1
14

MONITOR
1075
175
1200
220
Number of effectors
count turtles with [effectiveness = 1]
17
1
11

INPUTBOX
180
225
292
285
pop_limit
2000.0
1
0
Number

TEXTBOX
185
193
289
221
Maximum number\nof T cells:
11
0.0
1

INPUTBOX
209
629
320
689
CTLA4_time
1.0
1
0
Number

TEXTBOX
210
613
331
631
Effect on short contacts
11
0.0
1

INPUTBOX
208
543
270
603
K_prime
0.2
1
0
Number

INPUTBOX
287
543
347
603
K_time
0.2
1
0
Number

INPUTBOX
209
711
270
771
alfa
1.0
1
0
Number

INPUTBOX
207
474
272
534
mean_aff
0.02
1
0
Number

INPUTBOX
286
711
350
771
beta_2
0.02
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

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="max_div">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lambda">
      <value value="2.41E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-gradient?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Division_rad">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Prob_chem">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="what-a-boundary?">
      <value value="&quot;32 patches&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="what-a-output?">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thresh_level">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pop_limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_aff">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim_time">
      <value value="80320"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="life-time">
      <value value="7200"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
