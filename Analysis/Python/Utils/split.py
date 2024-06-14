import warnings
import numpy as np
from sklearn.model_selection._split import _BaseKFold, GroupsConsumerMixin
from sklearn.utils import check_random_state
from sklearn.utils.multiclass import type_of_target
from sklearn.utils.validation import check_array, column_or_1d

class GroupStratifiedKFold(GroupsConsumerMixin, _BaseKFold):
    """
    A modification of "StratifiedKFold" to perform KFold CV
    stratifying on a group variable. The modification is simply
    (i) adding the ability to use the existing group parameter and
    (ii) using the groups as stratification in place of y.
    """

    def __init__(self, n_splits=5, shuffle=False, random_state=None):
        super().__init__(
            n_splits=n_splits,
            shuffle=shuffle,
            random_state=random_state
        )

    def _make_test_folds(self, X, y, groups=None):
        rng = check_random_state(self.random_state)
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        groups = check_array(groups, input_name="groups", ensure_2d=False, dtype=None)
        type_of_target_groups = type_of_target(groups)
        allowed_target_types = ("binary", "multiclass")
        if type_of_target_groups not in allowed_target_types:
            raise ValueError(
                "Supported target types are: {}. Got {!r} instead.".format(
                    allowed_target_types, type_of_target_groups
                )
            )

        groups = column_or_1d(groups)

        _, groups_idx, groups_inv = np.unique(groups, return_index=True, return_inverse=True)
        _, class_perm = np.unique(groups_idx, return_inverse=True)
        groups_encoded = class_perm[groups_inv]

        n_classes = len(groups_idx)
        groups_counts = np.bincount(groups_encoded)
        min_groups = np.min(groups_counts)
        if self.n_splits > min_groups:
            warnings.warn(
                "The least populated group has only %d"
                " members, which is less than n_splits=%d."
                % (min_groups, self.n_splits),
                UserWarning,
            )

        groups_order = np.sort(groups_encoded)
        allocation = np.asarray(
            [
                np.bincount(groups_order[i :: self.n_splits], minlength=n_classes)
                for i in range(self.n_splits)
            ]
        )

        test_folds = np.empty(len(groups), dtype="i")
        for k in range(n_classes):
            folds_for_group = np.arange(self.n_splits).repeat(allocation[:, k])
            if self.shuffle:
                rng.shuffle(folds_for_group)
            test_folds[groups_encoded == k] = folds_for_group
        return test_folds

    def _iter_test_masks(self, X, y=None, groups=None):
        test_folds = self._make_test_folds(X, y, groups)
        for i in range(self.n_splits):
            yield test_folds == i