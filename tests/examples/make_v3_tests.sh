# To make the v3.0 files:

git clone https://github.com/MesserLab/SLiM.git && cd SLiM
git checkout v3.0
mkdir build_v3.0 && cd build_v3.0
cmake .. && make
cd ../..
SLiM/build_v3.0/slim recipe_nonWF.slim && mv recipe_nonWF.trees recipe_nonWF.v3.0.trees
SLiM/build_v3.0/slim recipe_WF.slim && mv recipe_WF.trees recipe_WF.v3.0.trees

# To make the v3.2 files:

git clone https://github.com/MesserLab/SLiM.git && cd SLiM
git checkout v3.2
mkdir build_v3.2 && cd build_v3.2
cmake .. && make
cd ../..
SLiM/build_v3.2/slim recipe_nonWF.slim && mv recipe_nonWF.trees recipe_nonWF.v3.2.trees
SLiM/build_v3.2/slim recipe_WF.slim && mv recipe_WF.trees recipe_WF.v3.2.trees
