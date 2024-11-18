#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2023-02-01
# Updates:
# =============================================================================
"""
Process data category
"""
import logging
logger_name = '.'.join(["PDB_EX", __name__])
logger = logging.getLogger(logger_name)


def checkKey(d_, l_key):
    """Check whether the keys are present in the dictionary.
    The values of keys don't matter, can be either of ('','.','?)

    Args:
        d_ (dict): dictionary to be checked
        l_key (list): keys to check

    Returns:
        bool: True if keys are present; False if not
    """
    if d_ and type(d_) == dict:
        for key in l_key:
            if key not in d_:
                return False
        else:
            return True
    else:
        return False

def checkElementsAllNonEmpty(l_):
    """Check whether an element is (' ','.','?)
    Args:
        l_ (list): list to check

    Returns:
        bool: True if all are non-empty
    """
    for each in l_:
        each = each.strip()
        if not each:
            return False
        elif each in ('.', '?'):
            return False
    return True

def mergeCatItemsUpdateValuesByRowIndex(d_to, d_from, tf_replace=False):
    """Merge/Combine two cat dicts only if the two cats have the SAME number of rows,
    taking the simple unions on items, and update values on each row based on row index only. 
    Resulting items/columns will be the simple union of the items/columns from both cats.
    Resulting values/rows will be based on d_to, with update from d_from
        if tf_replace is False:
            if the value in d_to is non-empty (not in ('','.','?')), it is taken to the resulting row;
            if the value in d_to is empty or non-existing, then the value from d_from is taken.
        if tf_replace is True:
            if the value is in d_from, use the value regarless of the value is empty or not

    e.g. d_to = {
        'item1':['value11t', 'value12t'], 
        'item2':['value21t', 'value22t']
        }
        * as data frame or matrix of:
        item1       item2
        value11t    value21t
        value12t    value22t

        d_from = {
        'item1':['value11f', 'value12f'],
        'item3':['value31f', 'value32f']
        }
        * as data frame or matrix of:
        item1       item3
        value11f    value31f
        value12f    value32f

    the merged union cat is:
    if tf_replace == True:
        d_combine = {
        'item1':['value11f', 'value12f'], 
        'item2':['value21t', 'value22t'],
        'item3':['value31f', 'value32f']
        }
        * as data frame or matrix of:
        item1       item2       item3
        value11f    value21t    value31f
        value12f    value22t    value32f
    if tf_replace == False:
        d_combine = {
        'item1':['value11t', 'value12t'], 
        'item2':['value21t', 'value22t'],
        'item3':['value31f', 'value32f']
        }
        * as data frame or matrix of:
        item1       item2       item3
        value11t    value21t    value31f
        value12t    value22t    value32f
        * If tf_replace is False, and value11t from d_to is in ('','.','?'), 
        then it's replaced by non-empty value11f from d_from

    Args:
        d_to (dict): primary dict, high priority
        d_from (dict): secondary dict to be merged to primary dict, low priority
        l_key2check (list): list of keys to check, with their values to be non-redundant for the combined cat

    Returns:
        dict: combined cat dict
    """
    if not d_from:
        return d_to
    if not d_to:
        return d_from

    # check whether both cats have the same number of rows, if not, return d_to
    if len(list(d_to.values())[0]) != len(list(d_from.values())[0]):
        return d_to

    # initialize combined cat
    l_key_to = list(d_to.keys())
    l_key_from = list(d_from.keys())
    l_key_combine = []
    for key_to in l_key_to:
        l_key_combine.append(key_to)
    for key_from in l_key_from:
        if key_from not in l_key_combine:
            l_key_combine.append(key_from)  # union the items from two cats
    d_combine = {}
    for key in l_key_combine:
        d_combine[key] = []

    # initialize d_combine with d_to values and union of items from two cats
    for i in range(len(list(d_to.values())[0])):
        for key in d_combine:
            if key in d_to:
                if d_to[key][i].strip():
                    d_combine[key].append(d_to[key][i])
                else:
                    d_combine[key].append('?')  # assign '?' for missing values
            else:
                d_combine[key].append('?')  # assign '?' for missing values

    # merge d_from to d_to by taking the union on values
    for i in range(len(list(d_from.values())[0])):
        for key in d_combine:
            if key in d_from:
                if tf_replace:
                    # update with value from d_from
                    d_combine[key][i] = d_from[key][i]
                else:
                    # if the value for an item in empty
                    if d_combine[key][i] in ('.', '?'):
                        if d_from[key][i].strip():
                            if d_from[key][i] not in ('.', '?'):
                                # update with non-empty value from d_from
                                d_combine[key][i] = d_from[key][i]

    return d_combine


def mergeCatItemsUpdateValuesByKeys(d_to, d_from, l_key2check, tf_replace=False):
    """Merge/Combine two cat dicts by taking the unions on items, and update values on rows with matching key values.
    Resulting items/columns will be the simple union of the items/columns from both cats.
    Resulting values/rows will be based on d_to, with update from d_from on rows with the matching unique tuple of values for l_key2check:
        if tf_replace is False:
            if the value in d_to is non-empty (not in ('','.','?')), it is taken to the resulting row;
            if the value in d_to is empty or non-existing, then the value from d_from is taken.
        if tf_replace is True:
            if the value is in d_from, use the value regarless of the value is empty or not

    e.g. d_to = {
        'key_item':['key_value1', 'key_value2'],
        'item1':['value11t, value12t'], 
        'item2':['value21t, value22t']
        }
        * as data frame or matrix of:
        key_item    item1       item2
        key_value1  value11t    value21t
        key_value2  value12t    value22t

        d_from = {
        'key_item':['key_value1', 'key_value3'],
        'item1':['value11f', 'value12f'],
        'item3':['value31f', 'value32f']
        }
        * as data frame or matrix of:
        key_item    item1       item3
        key_value1  value11f    value31f
        key_value3  value12f    value32f

    the merged union cat is:
        d_combine = {
        'key_item':[key_value1, key_value2],
        'item1':[value11t, value12t], 
        'item2':[value21t, value22t],
        'item3':[value31f, '?'     ]
        }
        * as data frame or matrix of:
        key_item    item1       item2       item3
        key_value1  value11t    value21t    value31f
        key_value2  value12t    value22t    ?
        * If tf_replace is False, and value11t from d_to is in ('','.','?'), 
        then it's replaced by non-empty value11f from d_from
        If tf_replace is True, then value11t will be replaced by value11f

    Args:
        d_to (dict): primary dict, high priority
        d_from (dict): secondary dict to be merged to primary dict, low priority
        l_key2check (list): list of keys to check, with their values to be non-redundant for the combined cat

    Returns:
        dict: combined cat dict
    """
    if not d_from:
        return d_to
    if not d_to:
        return d_from

    # check whether key items are present in both cats, return d_to if not
    if not checkKey(d_to, l_key2check):
        return d_to
    if not checkKey(d_from, l_key2check):
        return d_to

    # initialize combined cat
    l_key_to = list(d_to.keys())
    l_key_from = list(d_from.keys())
    l_key_combine = []
    for key_to in l_key_to:
        l_key_combine.append(key_to)
    for key_from in l_key_from:
        if key_from not in l_key_combine:
            l_key_combine.append(key_from)  # union the items from two cats
    d_combine = {}
    for key in l_key_combine:
        d_combine[key] = []

    # initialize d_combine with d_to values and union of items from two cats
    for i in range(len(list(d_to.values())[0])):
        for key in d_combine:
            if key in d_to:
                if d_to[key][i].strip():
                    d_combine[key].append(d_to[key][i])
                else:
                    d_combine[key].append('?')  # assign '?' for missing values
            else:
                d_combine[key].append('?')  # assign '?' for missing values

    # create dict of unique tuple values in d_combine (same as d_to after initialization)
    # with key as each unique tuple value and value as the row index
    d_value2check_to = {}
    for i in range(len(list(d_combine.values())[0])):
        l_row = []
        for key in l_key2check:
            l_row.append(d_combine[key][i].lower())
        t_row = tuple(l_row)
        d_value2check_to[t_row] = i

    # merge d_from to d_to by taking the union on values
    for i in range(len(list(d_from.values())[0])):
        # create key tuple value for each row to be checked
        l_row = []
        for key in l_key2check:
            l_row.append(d_from[key][i].lower())
        t_row = tuple(l_row)
        # check the tuple value against the unique tuple value list from d_to
        if t_row in d_value2check_to:  # merge rows from two cats with the same tuple value
            # get row index in d_combine
            i_row2update = d_value2check_to[t_row]
            for key in d_combine:
                if tf_replace:
                    if key in d_from:
                        # update with value from d_from
                        d_combine[key][i_row2update] = d_from[key][i]
                else:
                    # if the value for an item in empty
                    if d_combine[key][i_row2update] in ('.', '?'):
                        if key in d_from:
                            if d_from[key][i].strip():
                                if d_from[key][i] not in ('.', '?'):
                                    # update with non-empty value from d_from
                                    d_combine[key][i_row2update] = d_from[key][i]

    return d_combine


def mergeCatItemsUnionValuesByKeys(d_to, d_from, l_key2check, tf_replace=False):
    """Merge/Combine two cat dicts by taking the unions on both items and values.
    Resulting items/columns will be the simple union of the items/columns from both cats.
    Resulting values/rows will be the union based on the unique tuple of values for l_key2check:
        Rows with unique key tuple values from both cats will be separat rows (added from d_to first, then d_from):
            The missing values will be marked by '?'
        Rows with redundant key tuple values from both cats will be combined:
            For each item/column of a combined row, 
            if tf_replace is False:
                if the value in d_to is non-empty (not in ('','.','?')), it is taken to the resulting row;
                if the value in d_to is empty or non-existing, then the value from d_from is taken.
            if tf_replace is True:
                if the value is in d_from, use the value regarless of the value is empty or not

    e.g. d_to = {
        'key_item':['key_value1', 'key_value2'],
        'item1':['value11t', 'value12t'], 
        'item2':['value21t', 'value22t']
        }
        * as data frame or matrix of:
        key_item    item1       item2
        key_value1  value11t    value21t
        key_value2  value12t    value22t

        d_from = {
        'key_item':['key_value1', 'key_value3'],
        'item1':['value11f', 'value12f'],
        'item3':['value31f', 'value32f']
        }
        * as data frame or matrix of:
        key_item    item1       item3
        key_value1  value11f    value31f
        key_value3  value12f    value32f

    the merged union cat is:
        d_combine = {
        'key_item':['key_value1', 'key_value2', 'key_value3'],
        'item1':['value11t', 'value12t', 'value12f'], 
        'item2':['value21t', 'value22t', '?'     ],
        'item3':['value31f', '?',        'value32f']
        }
        * as data frame or matrix of:
        key_item    item1       item2       item3
        key_value1  value11t    value21t    value31f
        key_value2  value12t    value22t    ?
        key_value3  value12f    ?           value32f
        * If tf_replace is False, and value11t from d_to is in ('','.','?'), 
        then it's replaced by non-empty value11f from d_from
        If tf_replace is True, then value11t will be replaced by value11f

    Args:
        d_to (dict): primary dict, high priority
        d_from (dict): secondary dict to be merged to primary dict, low priority
        l_key2check (list): list of keys to check, with their values to be non-redundant for the combined cat

    Returns:
        dict: combined cat dict
    """
    if not d_from:
        return d_to
    if not d_to:
        return d_from

    # check whether key items are present in both cats, return d_to if not
    if not checkKey(d_to, l_key2check):
        return d_to
    if not checkKey(d_from, l_key2check):
        return d_to

    # initialize combined cat
    l_key_to = list(d_to.keys())
    l_key_from = list(d_from.keys())
    l_key_combine = []
    for key_to in l_key_to:
        l_key_combine.append(key_to)
    for key_from in l_key_from:
        if key_from not in l_key_combine:
            l_key_combine.append(key_from)  # union the items from two cats
    d_combine = {}
    for key in l_key_combine:
        d_combine[key] = []

    # initialize d_combine with d_to values and union of items from two cats
    for i in range(len(list(d_to.values())[0])):
        for key in d_combine:
            if key in d_to:
                if d_to[key][i].strip():
                    d_combine[key].append(d_to[key][i])
                else:
                    d_combine[key].append('?')  # assign '?' for missing values
            else:
                d_combine[key].append('?')  # assign '?' for missing values

    # create dict of unique tuple values in d_combine (same as d_to after initialization)
    # with key as each unique tuple value and value as the row index
    d_value2check_to = {}
    for i in range(len(list(d_combine.values())[0])):
        l_row = []
        for key in l_key2check:
            l_row.append(d_combine[key][i].lower())
        t_row = tuple(l_row)
        d_value2check_to[t_row] = i

    # merge d_from to d_to by taking the union on values
    for i in range(len(list(d_from.values())[0])):
        # create key tuple value for each row to be checked
        l_row = []
        for key in l_key2check:
            l_row.append(d_from[key][i].lower())
        t_row = tuple(l_row)
        # check the tuple value against the unique tuple value list from d_to
        if t_row in d_value2check_to:  # merge rows from two cats with the same tuple value
            # get row index in d_combine
            i_row2update = d_value2check_to[t_row]
            for key in d_combine:
                if tf_replace:
                    if key in d_from:
                        # update with value from d_from
                        d_combine[key][i_row2update] = d_from[key][i]
                else:
                    # if the value for an item in empty
                    if d_combine[key][i_row2update] in ('.', '?'):
                        if key in d_from:
                            if d_from[key][i].strip():
                                if d_from[key][i] not in ('.', '?'):
                                    # update with non-empty value from d_from
                                    d_combine[key][i_row2update] = d_from[key][i]
        else:  # add d_from row into d_combine if the tuple value is not in d_to
            l_keyvalue = []
            for key in l_key2check:
                l_keyvalue.append(d_from[key][i])
            if not checkElementsAllNonEmpty(l_keyvalue):
                continue
 
            for key in d_combine:
                if key in d_from:
                    d_combine[key].append(d_from[key][i])
                else:
                    d_combine[key].append('?')

    return d_combine
