# Falls_detection
This project aims to detect a fall of eldery who carry a sensor on their body in order to alert when a falling occurred. In this project I used machine learning tools in order to detect whether a record of 20 minutes of daily actions (such as walking around the house, driving, sitiing and typin ect.) contains a fall, and segments the minute that it occurred in. The records contains accelerometer, gyroscope ,pressure and temperature sensors.

functions:
main - the main code that uses the model on a new recored (new data) and detect the falls.
main_load_data - this code load the data and creates features for the main code.
Extract_features - this function gets the vectores of the acceleroeter, gyroscope, presure and temperature for every sample and creates features.
Retrain - this functions retrains the model on a new calssified dataset.
Retrain_load_data - loads the records for the Retarin function.
Calculate_CSF - the CSF calculation for feature selection.
Fixing_names - a function that fixes the names of the records in the folder.

veriables:
Model_RUSBOOST - my trained model.
Min_max_norm - the minimum and maximum values for the normalization of the feature matrix.
Best_features - the indices of the best features for the model.
Template_values - the values of the accelerometer template signal during falling, in order to check the correlation (feature).
