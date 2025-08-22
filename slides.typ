#import "@preview/mannot:0.3.0": *
#import "touying/lib.typ": *
#import "@preview/pinit:0.1.4": *
#import "@preview/xarrow:0.3.0": xarrow
#import "@preview/cetz:0.4.1"
#import "psi-slides-0.6.1.typ": *
#import "@preview/grayness:0.3.0": image-grayscale, image-transparency
#import "@preview/algorithmic:1.0.3"
#import algorithmic: algorithm

// color-scheme can be navy-red, blue-green, or pink-yellow
// #let s = psi-slides.register(aspect-ratio: "16-9", color-scheme: "pink-yellow")
#show: psi-theme.with(aspect-ratio: "16-9",
                      color-scheme: "pink-yellow",
                             config-info(
                                title: [Alternative linear response strategies],
                                subtitle: [Spin-resolution and orbital energy/ΔSCF equivalence],
                                author: [Edward Linscott],
                                date: datetime(year: 2025, month: 9, day: 25),
                                location: [Gandia],
                                references: [references.bib],
                             ))

#set footnote.entry(clearance: 0em)
#show bibliography: set text(0.6em)

#let primary = rgb("#dc005a")
#let secondary = rgb("#f0f500")

#let blcite(reference) = {
  text(fill: white, cite(reference))
}

#let delayedmark(start, content, tag: none, color: primary, mark: mark, color-before: black, alternatives: none) = {
   let entries = (mark(content, tag: tag, color: color-before),)*(start - 1) + (mark(content, tag: tag, color: color),)
   alternatives(repeat-last: true, ..entries)
}

#let methods-with-marks(self) = {
  let (uncover, only, alternatives) = utils.methods(self)
  let dm = delayedmark.with(alternatives: alternatives)
  (uncover, only, alternatives, dm, dm.with(mark: markhl, color-before: white))
}

#title-slide()

#matrix-slide(repeat: 3, self => [

  #let (uncover, only, alternatives) = utils.methods(self)

  #image("figures/figure_2_cropped_recoloured.svg", width: 100%)
spin-resolved linear response
#pause
],[
  #let data = read("figures/fig_en_curve_sl_annotated_zoom_recolored.svg", encoding: none)
  #alternatives()[
    #image-transparency(data, format: "svg", alpha: 100%)
    #text([orbital energy/ΔSCF equivalence \ (Koopmans functionals)])
  ][
    #image-transparency(data, format: "svg", alpha: 20%)
    #text(fill: gray, [orbital energy/ΔSCF equivalence \ (Koopmans functionals)])
  ]
])
= Spin-resolved linear response
==
#slide(repeat: 5, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

  #align(horizon,
    grid(columns: (1fr, 1fr), align: horizon + center, inset: 1em,
    [$ E_U = sum_(delayedmark(#2, I sigma, color: primary)) U^delayedmark(#3, I, color: secondary) / 2 "Tr"[bold(n)^(delayedmark(#2, I sigma, color: primary)) (1 - bold(n)^(delayedmark(#2, I sigma, color: primary)))] $],
    [$ U^delayedmark(#3, I, color: secondary) = [chi_0^(-1) - chi^(-1)]_(I I) $],
    [#pause #pause #pause functional treats spin channels separately],
    [#pause LR treats them together],
    )
  )

  #blcite(<Linscott2018>)
])

#focus-slide()[Why?]

==

#align(top,
  image("figures/prb2018_pdf.png", width: 100%, height: 100%)
)

== Spin-resolved linear response

#slide(repeat: 8, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

  #grid(columns: (1fr, auto, 1fr), align: horizon + center, inset: 1em,
  [conventional], [#sym.arrow.r], [spin-resolved],
[#pause $ chi_(I J) = (d n^I) / (d v^J)$],
[#pause #sym.arrow.r],
[$chi_(I J)^(sigma sigma') = (d n^(I sigma)) / (d v^(J sigma')) $],
[#pause $ mat(chi_(1 1), chi_(1 2); chi_(2 1), chi_(2 2))$],
[#pause #sym.arrow.r],
[$
mat(delayedmarkhl(#6, mat(delim: #none, chi^(arrow.t arrow.t)_(1 1), chi^(arrow.t arrow.b)_(1 1); chi^(arrow.b arrow.t)_(1 1), chi^(arrow.b arrow.b)_(1 1))), 
    mat(delim: #none, chi^(arrow.t arrow.t)_(1 2), chi^(arrow.t arrow.b)_(1 2); chi^(arrow.b arrow.t)_(1 2), chi^(arrow.b arrow.b)_(1 2));
    mat(delim: #none, chi^(arrow.t arrow.t)_(2 1), chi^(arrow.t arrow.b)_(2 1); chi^(arrow.b arrow.t)_(2 1), chi^(arrow.b arrow.b)_(2 1)), 
    mat(delim: #none, chi^(arrow.t arrow.t)_(2 2), chi^(arrow.t arrow.b)_(2 2); chi^(arrow.b arrow.t)_(2 2), chi^(arrow.b arrow.b)_(2 2)))
   
 $],
 [#pause #pause $U^I = [chi_0^(-1) - chi^(-1)]_(I I)$],
 [#pause #sym.arrow.r],
 [$U^(I sigma) = ???$],
  )


])

== What is screening $U$?

#grid(columns: (1fr, 1fr), align: horizon + center, inset: 0.5em,
cetz.canvas({
  import cetz.draw: *

  let counter = 0
  let positions = ((0, 0), (2, 5), (3, 2), (4, -1), (6, 2), (5, 5), (-1, 3),)
  for pos in positions {
    counter = counter + 1
    for i in range(1, counter) {
      line(name: "line", positions.at(i - 1), pos, stroke: gray + 1pt)
      // content("line", text(fill: gray, [$chi_(#i #counter)$]), midpoint: 0.5)
    }
  }
  counter = 0
  for pos in positions {
    counter = counter + 1
    circle(pos, radius: 0.5, fill: primary, stroke: none, name: "circle")
    content("circle", text(fill: white, [#counter]))
  }
}),
cetz.canvas({
  import cetz.draw: *

  rect((-4, -3), (4, 3), stroke: none, fill: secondary, alpha: 0.5)
  content((4, 3), [bath], anchor: "north-east", padding: 0.5em)
  circle((0, 0), radius: 0.5, fill: primary, stroke: none)
  
}),
  [all sites included in response matrix],
  [only one site included in response matrix],
  [bare $U$],
  [fully-screened $U$],
)

#slide(repeat: 11, self=> [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)
  #table(columns: (auto, auto, 1fr), align: horizon + center, inset: 1em,
    [
      #pause
      fully-screened],
    [
      #pause
      #set text(size: 0.8em)
      $
        mat(mat(delim: #none, chi^(arrow.t arrow.t)_(1 1),                            ;                            , chi^(arrow.b arrow.b)_(1 1)),;,
            mat(delim: #none, chi^(arrow.t arrow.t)_(2 2),                            ;                            , chi^(arrow.b arrow.b)_(2 2)))
      $ 
    ],
    [
      #pause
      $
        U^(I sigma) = 1 / (chi_0)^(sigma sigma)_(I I) - 1 / chi^(sigma sigma)_(I I)
      $
    ],

    [
      #pause
      screened by \ opposite spin],
    [
      #pause
      #set text(size: 0.8em)
      $
      mat(mat(delim: #none, chi^(arrow.t arrow.t)_(1 1), chi^(arrow.t arrow.b)_(1 1); chi^(arrow.b arrow.t)_(1 1), chi^(arrow.b arrow.b)_(1 1)),;,
          mat(delim: #none, chi^(arrow.t arrow.t)_(2 2), chi^(arrow.t arrow.b)_(2 2); chi^(arrow.b arrow.t)_(2 2), chi^(arrow.b arrow.b)_(2 2)))
      $ 
    ],
    [
      #pause
      $
      f^(sigma sigma')_I = [(chi_0)^(sigma sigma)_(I I)]^(-1) - [chi^(sigma sigma)_(I I)]^(-1) #pause \ f^(sigma sigma')_(I) stretch(arrow.r)^(???) U^I "or" U^(I sigma)
      $
    ],

    [
      #pause
      also screened by \ other Hubbard sites],
    [
      #pause
      #set text(size: 0.8em)
    $
    mat(mat(delim: #none, chi^(arrow.t arrow.t)_(1 1), chi^(arrow.t arrow.b)_(1 1); chi^(arrow.b arrow.t)_(1 1), chi^(arrow.b arrow.b)_(1 1)),
        mat(delim: #none, chi^(arrow.t arrow.t)_(1 2), chi^(arrow.t arrow.b)_(1 2); chi^(arrow.b arrow.t)_(1 2), chi^(arrow.b arrow.b)_(1 2));
        mat(delim: #none, chi^(arrow.t arrow.t)_(2 1), chi^(arrow.t arrow.b)_(2 1); chi^(arrow.b arrow.t)_(2 1), chi^(arrow.b arrow.b)_(2 1)), 
        mat(delim: #none, chi^(arrow.t arrow.t)_(2 2), chi^(arrow.t arrow.b)_(2 2); chi^(arrow.b arrow.t)_(2 2), chi^(arrow.b arrow.b)_(2 2)))
    $ 
    ],
    [
      #pause
      $ f^(sigma sigma')_(I J) = ... $ (left as an exercise to the reader)
    ]
  )
])

= Advantages of spin-resolved linear response
#focus-slide()[1. Conceptual consistency]
== 

#align(center + horizon, [spin-resolved linear response #sym.arrow.l.r spin-resolved DFT+_U_ functional])

#pause
#align(right + bottom, [... we didn't explore DFT+$U^sigma$; instead see BLOR@Burgess2023)])

#focus-slide()[2. Unconstrained constrained linear response]
== Unconstrained constrained linear response

Suppose we want to compute $ lr((d^2E_"Hxc") / (d (n^I)^2) |)_(mu^I)$

This is easy with spin-resolved LR:

$
  (d^2 E_"Hxc") / (d (n^I)^2) = & 1 / 2 (d v_"Hxc"^arrow.t + d v_"Hxc"^arrow.b) / d(n^arrow.t + n^arrow.b)
  = 1 / 2 (f^(arrow.t arrow.t ) d n^arrow.t + f^(arrow.t arrow.b) d n^arrow.b + f^(arrow.b arrow.t) d n^arrow.t + f^(arrow.b arrow.b) d n^arrow.b) / (d n^arrow.t + d n^arrow.b)
$

"Impose" the constraint by setting $d n^arrow.t = d n^arrow.b$ to get...

$
  lr((d^2E_"Hxc") / (d n^2) |)_(mu) = & 1/4 (f^(arrow.t arrow.t) + f^(arrow.b arrow.b) + f^(arrow.t arrow.b) + f^(arrow.b arrow.t))
$

This simple average is one choice (of many) for $M: f^(sigma sigma')_I arrow.r U^I$

#focus-slide()[3. We can recover conventional linear response]

== Conventional linear response

For conventional LR, $ v^(I arrow.t) = v^(I arrow.b)$ #pause, in which case:

$
d n = sum_sigma d n^(sigma) = sum_(sigma sigma') chi^(sigma sigma') d v^(sigma') = sum_(sigma sigma') chi^(sigma sigma') d v
arrow.double.r.long chi_"conv" = (d n) / (d v) = sum_(sigma sigma') chi^(sigma sigma')
$
#pause
Likewise,
$
  (epsilon^(-1))_"conv" = ... = 1 / 2 sum_(sigma sigma') (f chi)^(sigma sigma')
$
#pause
And thus

$
  U = (epsilon^(-1) - 1) chi^(-1) = 1 / 2 (sum_(sigma sigma') (f chi)^(sigma sigma')) / (sum_(sigma sigma') chi^(sigma sigma'))
$
#focus-slide()[4. $J$ is free]

==
As defined by 

$
  J = - 1 / 2 (d v_"Hxc"^arrow.t - d v_"Hxc"^arrow.b) / (d (n^arrow.t - n^arrow.b)) = - 1 / 4 ((f^(arrow.t arrow.t) - f^(arrow.b arrow.t)) d n^arrow.t - (f^(arrow.b arrow.b) - f^(arrow.t arrow.b)) d n^arrow.b)/ (d (n^arrow.t - n^arrow.b))
$

Different ways to define $J$:
+ while keeping $n = n^arrow.t + n^arrow.b$ fixed:
  $
    J = - 1 / 4 (f^(arrow.t arrow.t) - f^(arrow.b arrow.t) - f^(arrow.t arrow.b) + f^(arrow.b arrow.b))
  $
+ for a perturbation where $d v^arrow.t = - d v^arrow.b$

#pause
#focus-slide()[5. Easy to implement]

==
#slide(repeat: 5, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

  #align(center, [

  #alternatives()[
    #image("figures/fig_combined_U_J_and_occ.svg", height: 50%)
  ][
    #image("figures/fig_MnO_magmom.svg", height: 70%)
  ][
    #image("figures/fig_MnO_bandgap.svg", height: 70%)
  ]
  ])
  #v(-1em)
  #pause
  #pause
  #pause
#align(right, [... for more details see Linscott et al. 2018. #pause Since then used by many authors@Orhan2020@Lambert2023@MacEnulty2023@Moore2024@MacEnulty2024 and opened the door to DFT+$U$-inspired approaches@Burgess2023@Burgess2024])

])


#matrix-slide(repeat: 2, self => [

  #let (uncover, only, alternatives) = utils.methods(self)

  #let data = read("figures/figure_2_cropped_recoloured.svg", encoding: none)
  #alternatives()[
    #image-transparency(data, format: "svg", alpha: 100%)
    #text([spin-resolved linear response])
  ][
    #image-transparency(data, format: "svg", alpha: 20%)
    #text(fill: gray, [spin-resolved linear response])
  ]
],[
  #let data = read("figures/fig_en_curve_sl_annotated_zoom_recolored.svg", encoding: none)
  #alternatives()[
    #image-transparency(data, format: "svg", alpha: 20%)
    #text(fill: gray, [orbital energy/ΔSCF equivalence \ (Koopmans functionals)])
  ][
    #image-transparency(data, format: "svg", alpha: 100%)
    #text([orbital energy/ΔSCF equivalence \ (Koopmans functionals)])
  ]
])

= Koopmans functionals

== Total energy differences vs. eigenvalues

#align(horizon,
grid(align: horizon, columns: 2, column-gutter: 1em,
  [
We all know that DFT underestimates the band gap. But why? #pause

The exact Green's function has poles that correspond to total energy differences

$
  ε_i = cases(E(N) - E_i (N-1) & "if" i in "occ", E_i (N+1) - E(N) & "if" i in "emp")
$

#pause

but DFT does #emph[not]
],[
  #image(width: 20em, "figures/fig_en_curve_gradients_zoom_recolored.svg")
]
))

#focus-slide()[Core idea: impose this condition on DFT]

== Imposing generalised piecewise linearity
#align(horizon,
grid(align: horizon, columns: (1fr, auto), column-gutter: 1em,
  [Formally, every orbital $i$ should have an eigenenergy
  $
    epsilon_i^"Koopmans" = ⟨
      phi_i mid(|)hat(H)mid(|)phi_i
    ⟩ = frac(dif E, dif f_i)
  $
  that is
  - independent of $f_i$
  - equal to $Delta E$ of explicit electron addition/removal
],[
  #image(width: 20em, "figures/fig_en_curve_gradients_zoom_recolored.svg")
]
))

== Imposing generalised piecewise linearity
#slide(repeat: 3, self => [
  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

#align(horizon,
grid(align: horizon, columns: (1fr, auto), column-gutter: 1em,
[
  $
  E^"KI" &[rho, {rho_i}] =
  E^"DFT" [rho]
  \ & +
  sum_i (
    - delayedmarkhl(#2, integral_0^f_i lr(angle.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i angle.r) dif f, tag: #<remove_nonlin>, color: primary)
  \ &
    + delayedmarkhl(#3, f_i integral_0^1 lr(angle.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i angle.r) dif f, tag: #<restore_linear>, color: primary)
  )
$
// Bakes the total energy differences $E^"DFT" [rho^(f_i arrow.r 1)] - E^"DFT" [rho^(f_i arrow.r 0)]$ into the functional
#uncover("2-")[#annot(<remove_nonlin>, pos: bottom)[#align(center, [removes dependence on $f_i$])]]
#uncover("3-")[#annot(<restore_linear>, pos: bottom)[#align(center, [restores linear dependence on $f_i$])]]

],[
  #image(width: 20em, "figures/fig_en_curve_gradients_zoom_recolored.svg")
]))

])
== Comparison with DFT+_U_ (and BLOR)
#slide()[

#set table(
  fill: (x, y) =>
    if calc.rem(y, 2) == 0 { silver } else { silver.lighten(75%) },
    inset: 0.5em,
)
#show table.cell.where(x: 0): set text(style: "italic")

#show table.header: {
  set text(fill: primary, weight: "semibold")
}

#table(align: horizon, columns: (1fr, 2fr, 2fr), stroke: none,
    table.header([], [*DFT+_U_*], [*Koopmans*]),
    table.hline(),
   [seeks to correct...],
   uncover("2-")[erroneous curvature in total energies w.r.t. $N$],
   uncover("4-")[erroneous curvature in total energies w.r.t. $f_i forall i$],
   [in practice...],
   uncover("3-")[corrects curvature in total energies w.r.t. local manifold (BLOR does so more faithfully)],
   uncover("5-")[removes dependence of $epsilon_i$ on $f_i$ and guarantees $epsilon_i = E_i (N plus.minus 1) - E(N)$],
   [correction applied to...],
   [],
   [],
   [orbitals defined by...],
   [],
   [],
   [parametrised by...],
   [],
   [],
)
]
== Electronic screening via parameters
#slide(repeat: 3, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)
  $
    E^"KI" [{rho_i}] = &
    E^"DFT" [rho]
    +
    sum_i (
      - integral_0^f_i lr(angle.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i angle.r) dif f
      + f_i integral_0^1 lr(angle.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i angle.r) dif f
    )
    #pause
    \ = & E^"DFT" [rho]
    + sum_i {
      - (E^"DFT" [rho] - delayedmark(#3, E^"DFT" [rho^(f_i arrow.r 0)], tag: #<ENm1_hard>, color: primary))
      + f_i (delayedmark(#3, E^"DFT" [rho^(f_i arrow.r 1)], tag: #<ENp1_hard>, color: primary) - delayedmark(#3, E^"DFT" [rho^(f_i arrow.r 0)], tag: #<ENm1b_hard>, color: primary))
    }
    // uncover("5-", 
    // \ arrow.r E^"uKI" [{rho_i}] approx & 
    // E^"DFT" [rho]
    // \ & +
    // sum_i {
    //   - (E^"DFT" [rho] - delayedmark(#3, E^"DFT" [rho - rho_i], tag: #<ENm1>, color: primary))
    //   + f_i (delayedmark(#3, E^"DFT" [rho - rho_i + n_i], tag: #<ENp1>, color: primary) - delayedmark(#3, E^"DFT" [rho - rho_i], tag: #<ENm1b>, color: primary))
    // }
    // )
  $

  #pause
  #annot(<ENm1_hard>, pos: bottom)[#align(center, [cannot evaluate \ directly])]
  #annot(<ENp1_hard>, pos: bottom)[#align(center, [cannot evaluate \ directly])]
  #annot(<ENm1b_hard>, pos: bottom)[#align(center, [cannot evaluate \ directly])]
  #pause

  // Instead use a frozen-orbital picture:
  
  // $
  //  rho^(f_i arrow.r f)(bold(r)) approx rho(bold(r)) + (f - f_i) |phi^N_i (bold(r))|^2
  // $
  // 
  // very easy to evaluate -- but not at all accurate! Correct this _post hoc_ via a screening parameter i.e.
  // 
  // $
  //   E[rho^(f_i arrow.r f)] approx alpha_i E[rho + (f - f_i) |phi^N_i (bold(r))|^2]
  // $
])

#slide[
#align(center + horizon, 
  image("figures/fig_pwl_DFT.svg", height: 100%)
)
]
#slide[
#align(center + horizon, 
  image("figures/fig_pwl_uKI.svg", height: 100%)
)
]
#slide[
#align(center + horizon, 
  image("figures/fig_pwl_alphaKI.svg", height: 100%)
)
]

#slide(repeat: 5, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

$
  E^"KI"_bold(alpha) [rho, {rho_i}] approx & 
  E^"DFT" [rho]
  \ & +
  sum_i delayedmark(#3, alpha_i, tag: #<alpha>, color: primary) {
    - (E^"DFT" [rho] - delayedmark(#2, E^"DFT" [rho - rho_i], tag: #<ENm1>, color: primary))
    + f_i (delayedmark(#2, E^"DFT" [rho - rho_i + n_i], tag: #<ENp1>, color: primary) - delayedmark(#2, E^"DFT" [rho - rho_i], tag: #<ENm1b>, color: primary))
  }
$

#pause
#annot(<ENm1>, pos: bottom)[uses frozen orbitals]
#annot(<ENp1>, pos: bottom)[uses frozen orbitals]
#annot(<ENm1b>, pos: bottom)[uses frozen orbitals]
#pause
#annot(<alpha>, pos: bottom)[screening parameter]
#pause
which is easy to evaluate _e.g._
$ H^"KI"_(i j) = angle.l phi_j|hat(h)^"DFT" + alpha_i hat(v)_i^"KI"|phi_i angle.r #h(2cm) hat(v)^"KI"_i = - E_"Hxc" [rho - n_i] + E_"Hxc" [rho] - integral v_"Hxc" (bold(r)', [rho]) n_i d bold(r)' $

#pause
Screening parameters _not_ a fitting parameter!

])

// == Screening
// 
// #grid(columns: (1fr, 2fr), 
// [
// #align(center + horizon, 
//   image("figures/fig_pwl.png", width: 100%)
// )
// ],
// [
// 
// #pause
// Construct $alpha_i$ from explicit $Delta$SCF calculations@Nguyen2018@DeGennaro2022a
// 
// $
//   alpha_i = alpha_i^0 (Delta E_i - lambda_(i i)(0)) / (lambda_(i i)(alpha^0) - lambda_(i i)(0)) "where" lambda_(i i)(alpha) = angle.l phi_i|hat(h)^"DFT" + alpha hat(v)_i^"KI"|phi_i angle.r $
// 
// #pause
// Recast via linear response@Colonna2018:
// 
// $
//   alpha_i = (angle.l n_i mid(|) epsilon^(-1) f_"Hxc" mid(|) n_i angle.r) / (angle.l n_i mid(|) f_"Hxc" mid(|) n_i angle.r)
// $
// 
// which can be efficiently computed via DFPT@Colonna2022
// ],
// )

== Orbital-density dependence
#slide()[
The potential is orbital-density-dependent!
#v(-0.5em)
  $ v^"KI"_(i in"occ") = - E_"Hxc" [rho - n_i] + E_"Hxc" [rho] - integral v_"Hxc" (bold(r)', [rho]) n_i d bold(r)' $

#pause

- loss of unitary invariance@Nguyen2018
#v(-1em)
#align(center,
  grid(columns: (auto, auto), column-gutter: 1em,
  image("figures/fig_nguyen_variational_orbital.png", width: 10em),
  image("figures/fig_nguyen_canonical_orbital.png", width: 10em),
  [two variational orbitals],
  [a canonical orbital],
  )
) #pause
- we can use MLWFs@Marzari2012 #pause
- we know $hat(H)|phi_i angle.r$ but not $hat(H)$ #pause
- a natural generalisation of DFT towards spectral functional theory@Ferretti2014
]

== To summarise...
$
  E^"KI"_bold(alpha) [rho, {rho_i}] =
  E^"DFT" [rho] +
  sum_i alpha_i { &
    - (E^"DFT" [rho] - E^"DFT" [rho - rho_i])
  \ &
    + f_i (E^"DFT" [rho - rho_i + n_i] - E^"DFT" [rho - rho_i])
  }
$

- an orbital-by-orbital correction to DFT
- screening parameters
- orbital-density-dependence
- total energy at integer occupations unchanged!

== Comparison with DFT+_U_ (and BLOR)
#slide()[

#set table(
  fill: (x, y) =>
    if calc.rem(y, 2) == 0 { silver } else { silver.lighten(75%) },
    inset: 0.5em,
)

#show table.cell.where(y: 0): set text(weight: "bold", fill: primary)

#show table.cell.where(x: 0): set text(style: "italic")

#show table.cell: it => {
  set text(size: 0.9em)
  it
}

#table(align: horizon, columns: (1fr, 2fr, 2fr), stroke: none,
    table.header([], [DFT+_U_], [Koopmans]),
    table.hline(),
   [seeks to correct...],
   [erroneous global curvature in total energies w.r.t. $N$],
   [erroneous global curvature in total energies w.r.t. #uncover("6-")[*canonical*] orbital occupancies],
   [in practice...],
   [corrects curvature in total energies w.r.t. local manifold (BLOR does so more faithfully)],
   [removes dependence of $epsilon_i$ on #uncover("7-")[*variational*] orbital occupations and guarantees $epsilon_i = E_i (N plus.minus 1) - E(N)$],
   [correction applied to...],
   uncover("2-")[selected subspaces (e.g. _3d_ orbitals)],
   uncover("4-")[the entire system],
   [orbitals defined by...],
   uncover("3-")[Hubbard projectors (atom-centred, frozen, incomplete)],
   uncover("5-")[variational (localised) orbitals],
   [parametrised by...],
   uncover("8-")[${U_I}$], //, defined w.r.t. charge-neutral excitations if using LR],
   uncover("9-")[${alpha_i}$], // defined w.r.t. charged excitations]

)
  
]

= Results

== Molecular systems

=== Ionisation potentials@Colonna2019
#align(center + horizon,
image("figures/colonna_2019_gw100_ip.jpeg", width: 100%)
)

=== UV photoemission spectra@Nguyen2015
#align(center + horizon,
image("figures/fig_nguyen_prl_spectra_pink.png", width: 100%)
)


== Extended systems
#slide[
=== Prototypical semiconductors and insulators @Nguyen2018

#show table.cell: it => {
  if it.x == 3 or it.x == 4 {
    set text(fill: primary, weight: "semibold")
    it
  } else {
    it
  }
}

#grid(align: center + horizon, columns: 2, column-gutter: 1em,
image("figures/scatter_plot.png", height: 80%),
table(columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr), inset: 0.5em, stroke: none,
table.header([], [PBE], [G#sub[0]W#sub[0]], [KI], [KIPZ], [QSGW̃]),
table.hline(),
[$E_"gap"$], [2.54], [0.56], [0.27], [0.22], [0.18],
[IP], [1.09], [0.39], [0.19], [0.21], [0.49]
))
  
]

#slide[
=== ZnO @Colonna2022
#v(-1em)
#align(center + horizon,
grid(align: center + horizon, columns: 3, column-gutter: 1em,
image("figures/ZnO_lda_cropped.png", height: 50%),
image("figures/ZnO_hse_cropped_noaxis.png", height: 50%),
image("figures/ZnO_ki_cropped_noaxis.png", height: 50%),
))
#show table.cell: it => {
  set text(size: 0.8em)
  if it.x == 5 {
    set text(fill: primary, weight: "semibold")
    it
  } else {
    it
  }
}
#table(columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1.5fr), align: center, inset: 0.5em, stroke: none,
table.header([], [LDA ], [HSE ], [GW#sub[0] ], [scGW̃ ], [KI ], [exp ]),
table.hline(),
[$E_"gap"$], [0.79], [2.79], [3.0], [3.2], [3.68], [3.60],
[$angle.l epsilon_d angle.r$], [-5.1], [-6.1], [-6.4], [-6.7], [-6.93], [-7.5 to -8.81 ],
[$Delta$], [4.15], [], [], [], [4.99], [5.3]
)
  
]


= Caveats

== Limitations

- only valid for systems with $E_"gap"$ > 0 #pause
- empty state localisation in the bulk limit #pause
- can break crystal point group symmetry

== Resonance with other efforts

- Wannier transition state method of Anisimov and Kozhevnikov@Anisimov2005
- Optimally-tuned range-separated hybrid functionals of Kronik, Pasquarello, and others@Kronik2012@Wing2021
- Ensemble DFT of Kraisler and Kronik@Kraisler2013
- Koopmans-Wannier method of Wang and co-workers@Ma2016
- Dielectric-dependent hybrid functionals of Galli and co-workers@Skone2016a
- Scaling corrections of Yang and co-workers@Li2018


// == 
// #image("figures/supercell_workflow.png", width: 100%)
// 
// #image("figures/primitive_workflow.png", width: 65.5%)

#focus-slide()[
#align(center, image(width: 80%, "media/logos/koopmans_white_on_transparent.svg"))
]

== `koopmans`
#matrix-slide(alignment: horizon)[
  #image("figures/website_cropped.png")
][
  
  - automated workflows
  - `Quantum ESPRESSO` backend
  - easy installation
  - python API
  
  See `koopmans-functionals.org`
]

==
#align(center + horizon,
image("figures/supercell_workflow.png", width: 100%)
)

#matrix-slide(alignment: horizon, columns: (3fr, 2fr))[
  #image("figures/black_box_filled_square.png")
][
 
  Our goal:
  + accurate
  + robust
  + minimal input
  + fast

]


== 
#slide()[
#show table.cell: it => {
  if it.x == 5 {
    set text(fill: primary, weight: "semibold")
    it
  } else {
    it
  }
}
#grid(columns: 3, align: horizon + center,
  [
  #set text(size: 0.45em)
  #raw(read("scripts/gaas_auto.json"), block: true, lang: "json")
  ],
  [
    #set text(size: 3em)
    #sym.arrow.r
  ],
  [
    #set text(size: 0.7em)
    #image("figures/Unfold_And_Interpolate_bandstructure.png", height: 60%)
    #table(columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr), inset: 0.5em, stroke: none,
    table.header([], [LDA], [HSE], [GW#sub[0]], [scGW̃ ], [KI], [exp]),
    table.hline(),
    [$E_"gap"$], [0.26], [1.28], [1.55], [1.62], [1.54], [1.55],
    [$angle.l epsilon_d angle.r$], [-14.9], [-15.6], [-17.3], [-17.6], [-17.9], [-18.9],
    [$Delta$], [12.8], [13.9], [], [], [12.7], [13.1]
    )
  ]
)
]

= Summary
== Summary
#grid(
  columns: (1fr, 2fr),
  gutter: 1em,
  image("figures/black_box_filled_square.png", width: 100%),
  text[
    Koopmans functionals...
    - impose generalised piecewise linearity condition to DFT
    - give band structures with comparable accuracy to state-of-the-art GW
    - can be used in place of GW in BSE calculation of excitons, for systems with strong SOC, ...
    - are increasingly black-box
  ],
)

== Open questions

#pause
- why does correcting _local_ charged excitations correct the description of delocalized excitations? #pause
- is there a good metric for selecting variational orbitals (_i.e._ the subspace with respect to which we enforce piecewise linearity)? #pause
- are off-diagonal corrections appropriate? What form should they take? #pause
- how to extend to metallic systems? #pause
- can we provide a formal basis for the Koopmans correction? #pause
  - GKS
  - spectral functional theory@Ferretti2014
  - ensemble DFT
  - RDMFT

== Acknowledgements
#align(center + horizon, 
grid(columns: 9, column-gutter: 0.5em, align: center, row-gutter: 0.5em,
  image("media/mugshots/david_oregan.jpg", height: 40%),
  image("media/mugshots/andrew_burgess.jpeg", height: 40%),
  image("media/mugshots/nicola_colonna.png", height: 40%),
  image("media/mugshots/miki_bonacci.jpg", height: 40%),
  image("media/mugshots/aleksandr_poliukhin.jpg", height: 40%),
  image("media/mugshots/marija_stojkovic.jpg", height: 40%),
  image("media/mugshots/junfeng_qiao.jpeg", height: 40%),
  image("media/mugshots/yannick_schubert.jpg", height: 40%),
  image("media/mugshots/nicola_marzari.jpeg", height: 40%),
  [David O'Regan], [Andrew Burgess], [Nicola Colonna], [Miki Bonacci], [Aleksandr Poliukhin], [Marija Stojkovic], [Junfeng Qiao], [Yannick Schubert], [Nicola Marzari]
)
)

#align(
  center,
  grid(
    columns: 2,
    align: horizon + center,
    gutter: 2em,
    image("media/logos/SNF_logo_standard_web_color_pos_e.svg", height: 20%),
    image("media/logos/marvel_color_on_transparent.png", height: 20%),
  ),
)


#focus-slide()[#align(center, text(size: 2em, [Thank you!]) + linebreak() + text(size: 0.5em, style: "italic", [these slides are available at #h(0.2em) #box[#move(dy: 0.1em, image("media/logos/github-mark-white.svg", height: 1em))] `elinscott-talks`]))]


#show: appendix

#focus-slide()[#align(center, text(size: 2em, [spare slides]))]

== Frozen orbital approximation

  #v(-5em)
  #align(center + horizon, 
  grid(align: center + horizon, columns: 3, column-gutter: 2cm, row-gutter: 1cm,
  cetz.canvas({
    import cetz.draw: *
    content((1.25, 1.5), [$rho$])
    circle((0, 0), radius: 1, fill: primary, stroke: none)
    circle((2.5, 0), radius: 1, fill: primary, stroke: none)

  }),
  cetz.canvas({
    import cetz.draw: *

    content((9, 1.5), [$rho^(f_1 arrow.r 0)$])
    arc((10.75, 0), start: 0deg, stop: 360deg, radius: (1.5, 1), fill: primary, stroke: none)
    circle((8, 0), radius: 1, fill: none, stroke: (thickness: 2pt, paint: primary))
    circle((8, 0), radius: 1, fill: none, stroke: (dash: "dashed", thickness: 2pt, paint: white))
    // content((8, -1.5), [$f_1 = 0$])
  }),
  cetz.canvas({
    import cetz.draw: *

    content((17.25, 1.5), [$rho - |psi^N_1(r)|^2$])
    circle((16, 0), radius: 1, fill: none, stroke: (dash: "dashed", thickness: 2pt, paint: primary))
    circle((18.5, 0), radius: 1, fill: primary, stroke: none)
  }),
  [2-electron solution],
  [what we'd like to evaluate],
  [what we can quickly evaluate]

  ))

#matrix-slide(columns: (3fr, 2fr))[
#align(center + horizon,
  {only("1")[#image("figures/alpha_calc/fig_alpha_calc_step_0.png", height: 80%)]
  only("2")[#image("figures/alpha_calc/fig_alpha_calc_step_1.png", height: 80%)]
  only("3")[#image("figures/alpha_calc/fig_alpha_calc_step_2.png", height: 80%)]
  only("4-5")[#image("figures/alpha_calc/fig_alpha_calc_step_3.png", height: 80%)]
  only("6-7")[#image("figures/alpha_calc/fig_alpha_calc_step_4.png", height: 80%)]
  }
)
][
#only("7")[$ alpha_i = alpha_i^0 (Delta E_i - lambda_(i i)(0)) / (lambda_(i i)(alpha^0) - lambda_(i i)(0)) $
$ lambda_(i i)(alpha) = angle.l phi_i|hat(h)^"DFT" + alpha hat(v)_i^"KI"|phi_i angle.r $]
]

== Issues with extended systems

#align(center + horizon, 
  image("figures/fig_nguyen_scaling.png", width: 60%)
)

Two options: #pause _1._ use a more advanced functional#pause, or _2._ stay in the "safe" region
#blcite(<Nguyen2018>)

== 
#slide()[
#set text(size: 0.8em)
#raw(read("scripts/gaas.json"), block: true, lang: "json")
]
= Extensions
= Non-collinear spin
== Non-collinear spin

$ rho_i (bold(r)) pause arrow.r bold(rho)_i (bold(r)) = (rho_i (bold(r)), m_i^x (bold(r)), m_i^y (bold(r)), m_i^z (bold(r))) $

#pause

e.g. for the corrective potential

$ v_i^"qKI" = - 1 / 2 integral dif bold(r) dif bold(r)' rho_i (bold(r)) f_"Hxc" (bold(r), bold(r)') rho_i (bold(r)') + (1 - f_i) integral d bold(r)' f_"Hxc" (bold(r), bold(r)') rho_i (bold(r)') $

#pause

#align(center, sym.arrow.b)

$ v_i^"qKI" = - 1 / 2 integral dif bold(r) dif bold(r)' bold(rho)_i (bold(r)) bb(F)_"Hxc" (bold(r), bold(r)') bold(rho)_i (bold(r)') sigma_0 + (1 - f_i) sum_alpha integral d bold(r)' [bb(F)_"Hxc" (bold(r), bold(r)') bold(rho)_i (bold(r)')]_alpha sigma_alpha $


#blcite(<Marrazzo2024>)

#pagebreak()

CsPbBr#sub[3] #blcite(<Marrazzo2024>)
#v(-2em)
#align(center + horizon,
image("figures/marrazzo_CsPbBr3_bands.svg", height: 50%)
)
#table(align: center, columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1.5fr), inset: 0.5em, stroke: none,
table.header([], [LDA ], [HSE ], [G#sub[0]W#sub[0] ], [scGW̃ ], [*KI*], [exp ]),
table.hline(),
[*with SOC*], [0.18], [0.78], [0.94], [1.53], [*1.78*], [1.85],
[without SOC], [1.40], [2.09], [2.56], [3.15], [3.12], [],
)

= Optical spectra
== Optical spectra

Solve the BSE, using Koopmans eigenvalues in lieu of GW

#pause

#v(-1em)
#align(center + horizon,
grid(columns: 2,
image("figures/silicon_bse_spectra.png", height: 50%),
image("figures/silicon_bse_excitons.png", height: 50%)
))

#v(-1em)

#show table.cell: it => {
  set text(size: 0.8em)
  it
}
#table(align: center + horizon, columns: (auto, 1fr, 1fr, 1fr, 1fr), inset: 0.5em, stroke: none,
table.header([silicon], [indirect gap ], [direct gap ], [first excitonic peak ], [excitonic binding energy ]),
table.hline(),
[*qKI+BSE*], [1.12], [3.31], [3.42], [0.09], 
[G#sub[0]W#sub[0]+BSE], [1.17], [3.25], [3.34], [0.09],
)

= Computational cost and scaling
== Computational cost and scaling
#align(center + horizon,
image("figures/timings/benchmark.svg", width: 80%)
)

#pagebreak()

The vast majority of the computational cost: determining screening parameters

$
  alpha_i = (angle.l n_i|epsilon^(-1) f_"Hxc"|n_i angle.r) / (angle.l n_i|f_"Hxc"|n_i angle.r)
$

#pause

- a local measure of screening of electronic interactions #pause
- one screening parameter per orbital
- must be computed #emph[ab initio] via... #pause
  - $Delta$SCF@Nguyen2018@DeGennaro2022a: embarrassingly parallel steps which each cost $cal(O)(N_"SC"^3) tilde cal(O)(N_bold(k)^3 N^3)$ #pause
  - DFPT@Colonna2018@Colonna2022: $cal(O)(N_bold(k)^2 N^3)$

= Machine-learned electronic screening
== Machine-learned electronic screening

#pagebreak()

#slide[
  #align(
    center,
    grid(
      columns: 5,
      align: horizon,
      gutter: 1em,
      image("figures/orbital.emp.00191_cropped.png", height: 30%),
      $stretch(->)^("power spectrum decomposition")$,
      $vec(delim: "[", x_0, x_1, x_2, dots.v)$,
      $stretch(->)^("ridge regression")$,
      $alpha_i$,
    ),
  )

  $
    c^i_(n l m, k) & = integral dif bold(r) g_(n l) (r) Y_(l m)(theta,phi) n^i (
      bold(r) - bold(R)^i
    )
  $


  $
    p^i_(n_1 n_2 l,k_1 k_2) = pi sqrt(8 / (2l+1)) sum_m c_(n_1 l m,k_1)^(i *) c_(n_2 l m,k_2)^i
  $

  #blcite(<Schubert2024>)
]

#pagebreak()

#slide[
  #align(
    center,
    grid(
      columns: 2,
      align: horizon + center,
      gutter: 1em,
      image("figures/water.png", height: 70%),
      image("figures/CsSnI3_disordered.png", height: 70%),

      "water", "CsSnI" + sub("3"),
    ),
  )
  #blcite(<Schubert2024>)
]

The use-case

   #grid(columns: 8, column-gutter: 0.3em, row-gutter: 0.3em,
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        grid.cell(align: center + horizon, [...]),
        grid.cell(inset: 0.4em, align: center, fill: primary, colspan: 3, text(fill: white, "train", size: 1em, weight: "bold")),
        grid.cell(inset: 0.4em, align: center, fill: secondary, colspan: 5, text("predict", size: 1em, weight: "bold")),
  )

  #pause
  N.B. not a general model


#slide[
  #grid(
    columns: (1fr, 1fr),
    align: horizon + center,
    gutter: 1em,
    image(
      "figures/water_cls_calc_vs_pred_and_hist_bottom_panel_alphas.svg",
      height: 70%,
    ),
    image(
      "figures/CsSnI3_calc_vs_pred_and_hist_bottom_panel_alphas.svg",
      height: 70%,
    ),

    "water", "CsSnI" + sub("3"),
  )
  #blcite(<Schubert2024>)
]

#slide[
  #grid(
    columns: (1fr, 1fr),
    align: center + horizon,
    gutter: 1em,
    image(
      "figures/convergence_key.png",
      height: 5%,
    ) +  v(-1em) +
    image(
      "figures/convergence_fig.png",
      height: 55%,
    ),
    image("figures/speedup.png", height: 60%),

    [*accurate* to within $cal("O")$(10 meV) _cf._ typical band gap accuracy of $cal("O")$(100 meV)],
    [*speedup* of $cal("O")$(10) to $cal("O")$(100)],
  )

  #blcite(<Schubert2024>)
]

= Symmetries
== Taking advantage of symmetries
To compute screening parameters via DFPT...
#algorithm(inset: 0.3em, indent: 1em, {
  import algorithmic: *
  Function("CalculateAlpha", ($n$,), {
    For($bold(q) in "BZ"$,
    {
        For($bold(k) in "BZ"$, {Comment[Linear system $A x = b$ to obtain $Delta psi_(bold(k)+bold(q),v)(bold(r))$]})
          Assign[$Delta rho^(0n)_(q)$][$sum_(bold(k)v)psi^*_(bold(k)v) (bold(r))Delta psi_(bold(k)+bold(q),v)(bold(r)) + c.c.$]
          Assign[$Pi^((r))_(0 n, bold(q))$][$angle.l Delta rho^(0 n)_(bold(q))|f_"Hxc"|rho^(0 n)_(bold(q)) angle.r$]
          Assign[$Pi^((u))_(0 n, bold(q))$][$angle.l rho^(0 n)_bold(q)|f_"Hxc"|rho^(0 n)_bold(q) angle.r$]
    })
    Return[$1 + sum_bold(q) Pi^((r))_(0 n, bold(q)) \/ sum_bold(q) Pi^((u))_(0 n, bold(q))$]
  })
})

#pagebreak()

#align(center,
  image("figures/bz-to-ibz-outer.svg", height: 80%)
)
$bold(q) in "BZ" $ $arrow.r$ $bold(q) in "IBZ"(n)$ (the symmetry of the perturbation; lower than that of the primitive cell)
#pagebreak()
#align(center,
  image("figures/bz-to-ibz-inner.svg", height: 80%)
)
$bold(k) in "BZ"$ $arrow.r$ $bold(k) in "IBZ"(bold(q))$ (can only use symmetries that leave $bold(q)$ invariant)

#align(horizon + center, image("figures/bz-to-ibz-speedup.svg", height: 100%))

= Automated Wannierisation
== Automated Wannierisation
#slide()[
  Koopmans functionals rely heavily on Wannier functions...
  - to initialise the minmising orbitals, _or_
  - in place of the minimising orbitals entirely

#pause

#grid(
  columns: (2fr, 2fr, 3fr),
  align: center + horizon,
  gutter: 1em,
  image("figures/proj_disentanglement_fig1a.png", height: 45%),
  image("figures/new_projs.png", height: 45%),
  image("figures/target_manifolds_fig1b.png", height: 45%),

  text("projectability-based disentanglement") + cite(<Qiao2023>),
  text("use PAOs found in pseudopotentials"),
  text("parallel transport to separate manifolds") + cite(<Qiao2023a>),
)
]

== 
#blcite(<Huber2020>)
#v(-2em)
#align(center,
  [
  #grid(columns: 3, align: horizon, column-gutter: 0.5em,
    image("media/logos/koopmans_grey_on_transparent.svg", height: 3em),
    image("figures/handshake.png", height: 2em, alt: "handshake"),
    image("media/logos/aiida.svg", height: 3em)
  )
  #pause `$ koopmans run tio2.json` #pause $arrow.r$ `$ koopmans run --engine=aiida tio2.json`
  ]
)

remote compute, parallel step execution, provenance-tracking, (requires configuration, WIP...)

#pause
#align(center, 
  image("figures/aiida-speed-up.svg", width: 70%)
)


== Connections with approx. self-energies

#blcite(<Ferretti2014>)#blcite(<Colonna2019>)

Orbital-density functional theory:

$ (h + alpha_i v^(K I)_i)|psi_i angle.r = lambda_i|psi_i angle.r $ $v_i^(K I)(bold(r))$ is real, local, and state-dependent #pause

cf. Green's function theory:

$ (h + Sigma_i)|psi_i angle.r = z_i|psi_i angle.r $ $Sigma_i (bold(r), bold(r)')$ is complex, non-local, and state-dependent

#slide[
Hartree-Fock self-energy in localized representation

$Sigma_x (bold(r), bold(r)') = - sum_(k sigma)^("occ") psi_(k sigma)(bold(r)) & f_H (bold(r), bold(r'))psi^*_(k sigma)(bold(r)') \
& arrow.r.double.long angle.l phi_(i sigma)|Sigma_x|phi_(j sigma') angle.r approx - angle.l phi_(i sigma)|v_H [n_(i sigma)]|phi_(i sigma)angle.r delta_(i j)delta_(sigma sigma')$

Unscreened KIPZ#sym.at Hartree ($v_"xc" arrow.r 0$; $f_"Hxc" arrow.r f_H$; $epsilon^(-1) arrow.r 1$)

$angle.l phi_(i sigma)|v^"KIPZ"_(j sigma',"xc")|phi_(j sigma') angle.r
approx {(1/2 - f_(i sigma)) angle.l n_(i sigma)|f_H|n_(i sigma) angle.r - E_H [n_(i sigma)]}
approx - angle.l phi_(i sigma)|v_H [n_(i sigma)]|phi_(i sigma)angle.r delta_(i j)delta_(sigma sigma')$

]

#slide[
Screened exchange plus Coulomb hole (COHSEX)

$ Sigma^"SEX"_"xc" (bold(s), bold(s)') = - sum_(k sigma)^"occ" psi_(k sigma)(bold(r)) psi_(k sigma)^*(bold(r)) W(bold(r), bold(r)') $

$ Sigma^"COH"_"xc" (bold(s), bold(s)') = 1/2 delta(bold(s), bold(s)'){W(bold(r), bold(r)') - f_H (bold(r), bold(r)')} $

$ arrow.r.double.long angle.l phi_(i sigma)|Sigma^"COHSEX"_"xc"|phi_(j sigma')angle.r approx {(1/2 - f_(i sigma)) angle.l n_(i sigma)|W|n_(i sigma)angle.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$

KIPZ#sym.at Hartree with RPA screening ($v_"xc" arrow.r 0$; $f_"Hxc" arrow.r f_H$; $epsilon^(-1) arrow.r "RPA"$)

$ angle.l phi_(i sigma)|v^"KIPZ"_(j sigma',"xc")|phi_(j sigma')angle.r approx{(1/2 - f_(i sigma)) angle.l n_(i sigma)|W|n_(i sigma)angle.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$
]

#slide[
  Static GWΓ#sub[xc] --- local (DFT-based) vertex corrections@Hybertsen1987@DelSole1994

  $ Sigma^(G W Gamma_"xc")_"xc"(1, 2) = i G(1, 2) W_(t-e) (1, 2) $
  
  $ W_(t-e) = (1 - f_"Hxc" chi_0)^(-1) f_H $

  $ arrow.r.double.long angle.l phi_(i sigma)|Sigma^(G W Gamma_"xc")_"xc"|phi_(j sigma')angle.r approx{(1/2 - f_(i sigma)) angle.l n_(i sigma)|W_(t-e)|n_(i sigma)angle.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$

  KIPZ#sym.at DFT ($v_"xc" arrow.r$ DFT; $f_"Hxc" arrow.r$ DFT; $epsilon^(-1) arrow.r$ DFT)

  $ angle.l phi_(i sigma)|v^"KIPZ"_(j sigma',"xc")|phi_(j sigma')angle.r approx{angle.l phi_(i sigma)|v^"DFT"_(sigma,"xc")|phi_(i sigma)angle.r + (1/2 - f_(i sigma)) angle.l n_(i sigma)|epsilon^(-1)_(t-e) f_"Hxc"|n_(i sigma)angle.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$
]

== Model systems
=== Hooke's atom@Schubert2023

#align(center + horizon, 
  image("figures/schubert_vxc_only.jpeg", height: 70%)
)

= References
== References
#slide()[
#show bibliography: set text(0.95em)
#bibliography("references.bib", style: "nature-footnote.csl", title: none)

]