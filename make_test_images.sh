#!/usr/bin/fish

set destination_dir sample_images
set xdim 64
set ydim 64
set locations \
-0.04524074130409,0.9868162207157838,4.6E-12 \
-0.0452407411,0.9868162204352258,5.4E-10
set location_names swirls mandelflower

set magma_color_map \
"#FCFFB2" \
"#FCDF96" \
"#FBC17D" \
"#FBA368" \
"#FA8657" \
"#F66B4D" \
"#ED504A" \
"#E03B50" \
"#C92D59" \
"#B02363" \
"#981D69" \
"#81176D" \
"#6B116F" \
"#57096E" \
"#43006A" \
"#300060" \
"#1E0848" \
"#110B2D" \
"#080616" \
"#000005"

set magma_colortable $destination_dir/magma_colortable.png
convert xc:{$magma_color_map} +append $magma_colortable

set gamma 1
set options -c -x $xdim -y $ydim

for i in (seq (count $location_names))
  set output_filename $destination_dir/$location_names[$i]
  ./mandelbars $options -L $locations[$i] $output_filename.pgm
  convert $output_filename.pgm -gamma $gamma $magma_colortable -clut $output_filename.png
  rm $output_filename.pgm
end

rm $magma_colortable
