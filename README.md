# cnn-hep-thesis
Undergrad thesis. Using  a CNN to classify different HEP processes. 
For this repository to be useful we recommend to previously have some experience with CERN Open data. You can find their last workshop [here](https://cms-opendata-workshop.github.io/2022-08-01-cms-open-data-workshop/). \ 

## src 

The code you can found in this section (src) is meant to work in three different programs. The analyzers and python configuration files run in a cmssw docker environment. The image_generator_imass_script works in any notebook with the appropiate dependencies and finally the cnnmodel file works with google collab notebooks but it should work locally too. Due to the lack of gpu resources we trained the model using the google services, but this could be useful for any person trying to replicate the results given that the datasets are publicly avaible. 

