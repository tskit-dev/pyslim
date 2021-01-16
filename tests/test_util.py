"""
Test cases for utility functions.
"""
import tests
import numpy as np

import pyslim

class TestUniqueLabelsByGroup():

    def verify_unique_labels_by_group(self, group, label, minlength):
        x = pyslim.util.unique_labels_by_group(group, label, minlength)
        assert len(x) >= minlength
        for g in range(len(x)):
            u = set(label[group == g])
            if (len(u) == 1) != x[g]:
                print(g, u, x[g], label[group == g], label.dtype)
            assert (len(u) == 1) == x[g]

    def test_all_same(self):
        n = 10
        group = np.repeat(1, 10)
        label = np.arange(10)
        self.verify_unique_labels_by_group(group, label, 1)
        x = pyslim.util.unique_labels_by_group(group, label, 1)
        assert len(x) == 2
        assert x[0] == False
        assert x[1] == False
        label = np.repeat(5, 10)
        self.verify_unique_labels_by_group(group, label, 1)
        x = pyslim.util.unique_labels_by_group(group, label, 1)
        assert len(x) == 2
        assert x[0] == False
        assert x[1] == True

    def test_all_unique(self):
        ng = 10
        group = np.arange(ng)
        label = np.arange(ng)
        self.verify_unique_labels_by_group(group, label, ng)
        x = pyslim.util.unique_labels_by_group(group, label, ng)
        assert np.all(x)
        group = np.append(group, [-1, -1, -1])
        label = np.append(label, [0, 1, 2])
        self.verify_unique_labels_by_group(group, label, ng)
        x = pyslim.util.unique_labels_by_group(group, label, ng)
        assert np.all(x)

    def test_unique_labels_by_group(self):
        np.random.seed(23)
        for ng in 3 * np.arange(2, 15):
            for n in (10, 100):
                for nl in (2, ng):
                    for minl in (-5, 10000000):
                        group = np.random.choice(np.arange(ng) - 1, size=n)
                        # "integer" labels
                        label = minl + np.random.choice(np.arange(nl), size=n)
                        self.verify_unique_labels_by_group(group, label, ng)
                        # int32 labels
                        self.verify_unique_labels_by_group(group, label.astype("int32"), ng)
                        # and float labels
                        label = minl + np.random.choice(np.random.uniform(0, 1, nl), size=n)
                        self.verify_unique_labels_by_group(group, label, ng)


