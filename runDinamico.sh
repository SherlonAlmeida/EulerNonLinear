echo "Compilando programa..."
make

echo " "
echo "Gerando o arquivo GNUPLOT..."
python3 geraPlotagem.py

clear

echo " "
echo "Executando o Calculo de Combustao..."
./euler

echo " "
echo "Executando o GPROF..."
gprof euler &> profilingGPROF.dat

echo " "
echo "Gerando as imagens..."
gnuplot vorticidadeDinamico.gnu

echo " "
echo "Gerando o GIF..."
convert -delay 20 *.jpeg -loop 0 saida_animada.gif
if [ ! -d "Resultados" ]; then
    mkdir Resultados
else
    rm -r Resultados
    mkdir Resultados
fi
mv *.jpeg saida_animada.gif *.dat Resultados

make clean
