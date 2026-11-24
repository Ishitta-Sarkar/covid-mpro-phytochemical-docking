reinitialize
bg_color white
set ray_opaque_background, 1
load structures/receptor/6LU7_chainA.pdb, receptor
load structures/ligands/n3_crystal.pdb, crystal
load structures/poses/n3_top.pdb, redocked
hide everything
hide (elem H)
show cartoon, receptor
show sticks, crystal or redocked
set stick_radius, 0.18
color grey70, receptor
color cyan, crystal
color salmon, redocked
select site, chain A and (resi 41 or resi 145)
show sticks, site
color yellow, site
orient crystal
zoom crystal, 8
