import json
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model, datasets
from sklearn.metrics import mean_squared_error, r2_score


def init():
    global logreg
    config = json.load(open("config.json"))
    train_data = json.load(open("train_data.json"))
    test_data = json.load(open("test_data.json"))

    X_train = []
    X_test = []
    Y_train = []
    Y_test = []

    for set in train_data["sets"]:
        X_train.append(normal_data(set, config))
        Y_train.append(int(set["score"]))

    for set in test_data["sets"]:
        X_test.append(normal_data(set, config))
        Y_test.append(int(set["score"]))

    h = .02  # step size in the mesh
    print(len(X_train))
    logreg = linear_model.LogisticRegression(C=1e5)

    # we create an instance of Neighbours Classifier and fit the data.
    logreg.fit(X_train, Y_train)
    Y_pred = logreg.predict(X_test)
    print(X_test)
    Y_pred_proba = logreg.predict_proba(X_test)
    #print(int(Y_pred_proba[][1]*100))
    Y_pred = list(Y_pred)
    #print("Real val:",Y_test)
    #print("Pred val:",Y_pred)
    print("Score logistic regression score:")
    print("Mean squared error: %.2f" % mean_squared_error(Y_test, Y_pred))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.2f' % r2_score(Y_test, Y_pred))


def normal_data(set, config):
    data_vector = []

    f_tm = (float(set["F_tm"])-float(config["Tm"]["Min"]))/(float(config["Tm"]["Max"])-float(config["Tm"]["Min"]))
    data_vector.append(f_tm)
    r_tm = (float(set["R_tm"])-float(config["Tm"]["Min"]))/(float(config["Tm"]["Max"])-float(config["Tm"]["Min"]))
    data_vector.append(r_tm)

    f_gc = (float(set["F_gc"]) - float(config["GC Percent"]["Min"])) / (float(config["GC Percent"]["Max"]) - float(config["GC Percent"]["Min"]))
    data_vector.append(f_gc)
    r_gc = (float(set["R_gc"]) - float(config["GC Percent"]["Min"])) / (float(config["GC Percent"]["Max"]) - float(config["GC Percent"]["Min"]))
    data_vector.append(r_gc)

    f_len = (float(set["F_len"]) - float(config["Length"]["Min"])) / (float(config["Length"]["Max"]) - float(config["Length"]["Min"]))
    data_vector.append(f_len)
    r_len = (float(set["R_len"]) - float(config["Length"]["Min"])) / (float(config["Length"]["Max"]) - float(config["Length"]["Min"]))
    data_vector.append(r_len)

    tm_dif = (float(set["Tm_dif"]) - float(config["Temperature difference"]["Min"])) / (float(config["Temperature difference"]["Max"]) - float(config["Temperature difference"]["Min"]))
    data_vector.append(tm_dif)

    amp_len = (float(set["Amp_len"]) - float(config["Amplicon Length"]["Min"])) / (float(config["Amplicon Length"]["Max"]) - float(config["Amplicon Length"]["Min"]))

    data_vector.append(amp_len)

    return data_vector

def normal_set(set, config):
    data_vector = []

    f_tm = (float(set.forward_primer.primer_tm())-float(config["Tm"]["Min"]))/(float(config["Tm"]["Max"])-float(config["Tm"]["Min"]))
    data_vector.append(f_tm)
    r_tm = (float(set.reverse_primer.primer_tm())-float(config["Tm"]["Min"]))/(float(config["Tm"]["Max"])-float(config["Tm"]["Min"]))
    data_vector.append(r_tm)

    f_gc = (float(set.forward_primer.precent_gc()) - float(config["GC Percent"]["Min"])) / (float(config["GC Percent"]["Max"]) - float(config["GC Percent"]["Min"]))
    data_vector.append(f_gc)
    r_gc = (float(set.reverse_primer.precent_gc()) - float(config["GC Percent"]["Min"])) / (float(config["GC Percent"]["Max"]) - float(config["GC Percent"]["Min"]))
    data_vector.append(r_gc)

    f_len = (float(set.forward_primer.length) - float(config["Length"]["Min"])) / (float(config["Length"]["Max"]) - float(config["Length"]["Min"]))
    data_vector.append(f_len)
    r_len = (float(set.reverse_primer.length) - float(config["Length"]["Min"])) / (float(config["Length"]["Max"]) - float(config["Length"]["Min"]))
    data_vector.append(r_len)

    tm_dif = (float(set.tm_dif()) - float(config["Temperature difference"]["Min"])) / (float(config["Temperature difference"]["Max"]) - float(config["Temperature difference"]["Min"]))
    data_vector.append(tm_dif)

    amp_len = (float(set.get_amplicon_length()) - float(config["Amplicon Length"]["Min"])) / (float(config["Amplicon Length"]["Max"]) - float(config["Amplicon Length"]["Min"]))

    data_vector.append(amp_len)

    return data_vector

