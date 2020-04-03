mkdir tmp_data
K=50
for i in {1..50};
	do
	cp -r /home/tobias/Documents/Work/CompBio/PR/debug-data/HybridLinker_2015pathlinker-weighted_Wnt_10000 tmp_data/"$i"_2015pathlinker-weighted_Wnt_10000;
	python3 /home/tobias/Documents/Work/CompBio/PR/make_pr.py tmp_data "$i"_2015pathlinker-weighted_Wnt_10000;
done

for i in {1..50};
	do
	cp tmp_data/1_2015pathlinker-weighted_Wnt_10000/negatives.csv tmp_data/"$i"_2015pathlinker-weighted_Wnt_10000/negatives.csv;
done

mkdir dummy_plot
python3 /home/tobias/Documents/Work/CompBio/PR/plot_pr.py tmp_data dummy_plot $(ls tmp_data)

