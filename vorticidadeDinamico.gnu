reset
set term jpeg enhanced large font times 20 lw 2 size 1200,800
#set term jpeg enhanced  size 1200,1200
#set key at 2.4,1.35
#set key at 2.4,1.35
set key samplen 1
set key spacing 1
#set key font ",20"
set key outside
set zero 1e-12
set size ratio 0.2
set palette defined ( 0 "#7f0000",\
                      1 "#ee0000",\
                      2 "#ff7000",\
                      3 "#ffee00",\
                      4 "#90ff70",\
                      5 "#0fffee",\
                      6 "#0090ff",\
                      7 "#000fff",\
                      8 "#000090")
set cbrange [-2.0:1.0]
#set cbrange [0.5:4.2]
#set cbrange [-10:0]
#set xlabel font '100'
#set cblabel font '100'
#show cblabel
#set contour base
set cntrparam levels auto 20
# discrete 0.06,0.05,0.04,0.03,0.02,0.01,0.001,-0.01,-0.02,-0.03,-0.04,-0.05 #
#set cntrparam levels  increment 0.80,0.005,1.05
#set cntrparam levels  increment 0.80,0.005,1.05
#set cbtics 0.05
#set xrange[044:48]
set xrange[00:15]
set yrange[-2.5:2.5]
set isosamples 10000
set hidden3d
set pm3d
#set dgrid3d
set view map
unset surface
set nokey
set title 'Vorticity Countours Subsonic Mixing Layer'
#set key horizontal outside
set output 'subsonic0.jpeg'
splot 'contorno_0.dat'  u 1:2:7  w pm3d
