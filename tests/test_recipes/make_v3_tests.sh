# To make example files for previous file versions

git clone https://github.com/MesserLab/SLiM.git && cd SLiM

# To make the v3.0 files:

git checkout v3.0
mkdir -p build_v3.0 && cd build_v3.0
cmake .. && make
cd ../..
SLiM/build_v3.0/slim recipe_nonWF.slim && mv recipe_nonWF.trees recipe_nonWF.v3.0.trees
SLiM/build_v3.0/slim recipe_WF.slim && mv recipe_WF.trees recipe_WF.v3.0.trees

# To make the v3.2 files:

git checkout v3.2
mkdir -p build_v3.2 && cd build_v3.2
cmake .. && make
cd ../..
SLiM/build_v3.2/slim recipe_nonWF.slim && mv recipe_nonWF.trees recipe_nonWF.v3.2.trees
SLiM/build_v3.2/slim recipe_WF.slim && mv recipe_WF.trees recipe_WF.v3.2.trees

# To make the v3.3.1 files:

git checkout v3.3.1
mkdir -p build_v3.3.1 && cd build_v3.3.1
cmake .. && make
cd ../..
SLiM/build_v3.3.1/slim recipe_nonWF.slim && mv recipe_nonWF.trees recipe_nonWF.v3.3.1.trees
SLiM/build_v3.3.1/slim recipe_WF.slim && mv recipe_WF.trees recipe_WF.v3.3.1.trees

# To make the v3.4 files:

TAG=v3.4
git checkout $TAG
mkdir -p build_$TAG && cd build_$TAG
cmake .. && make
cd ../..
SLIMDIR=SLiM/build_${TAG}
$SLIMDIR/slim recipe_nonWF.slim && mv recipe_nonWF.trees recipe_nonWF.${TAG}.trees
$SLIMDIR/slim recipe_WF.slim && mv recipe_WF.trees recipe_WF.${TAG}.trees
git add -f recipe_nonWF.${TAG}.trees recipe_WF.${TAG}.trees

# To make the v3.5 files:

TAG=v3.5
git checkout $TAG
mkdir -p build_$TAG && cd build_$TAG
cmake .. && make
SLIMDIR=$(pwd)
cd ../..
$SLIMDIR/slim recipe_nonWF.slim && mv out.trees recipe_nonWF.${TAG}.trees
$SLIMDIR/slim recipe_WF.slim && mv out.trees recipe_WF.${TAG}.trees
git add -f recipe_nonWF.${TAG}.trees recipe_WF.${TAG}.trees

# To make the v3.6 files:

TAG=v3.6
git checkout $TAG
mkdir -p build_$TAG && cd build_$TAG
cmake .. && make
SLIMDIR=$(pwd)
cd ../..
$SLIMDIR/slim recipe_nonWF.slim && mv out.trees recipe_nonWF.${TAG}.trees
$SLIMDIR/slim recipe_WF.slim && mv out.trees recipe_WF.${TAG}.trees
git add -f recipe_nonWF.${TAG}.trees recipe_WF.${TAG}.trees

# To make the "mixed" files:

python3 make_v3_test_additions.py

# To make the v3.7 files:

TAG=v3.7
git checkout $TAG
mkdir -p build_$TAG && cd build_$TAG
cmake .. && make
SLIMDIR=$(pwd)
cd ../..
$SLIMDIR/slim recipe_nonWF.slim && mv out.trees recipe_nonWF.${TAG}.trees
$SLIMDIR/slim recipe_WF.slim && mv out.trees recipe_WF.${TAG}.trees
git add -f recipe_nonWF.${TAG}.trees recipe_WF.${TAG}.trees

# To make the v4.2.2 files:

TAG=v4.2.2
git checkout $TAG
mkdir -p build_$TAG && cd build_$TAG
cmake .. && make
SLIMDIR=$(pwd)
cd ../..
$SLIMDIR/slim recipe_nonWF.slim && mv out.trees recipe_nonWF.${TAG}.trees
$SLIMDIR/slim recipe_WF.slim && mv out.trees recipe_WF.${TAG}.trees
$SLIMDIR/slim recipe_WF_X.slim && mv out.trees recipe_WF_X.${TAG}.trees
$SLIMDIR/slim recipe_WF_Y.slim && mv out.trees recipe_WF_Y.${TAG}.trees
git add -f recipe_nonWF.${TAG}.trees recipe_WF.${TAG}.trees recipe_WF_X.${TAG}.trees recipe_WF_Y.${TAG}.trees


