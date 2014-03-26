#!/usr/local/uvcdat/latest/bin/python
#
# =========================
#  EG python support procs
# =========================
#

def find_element_in_list(element, list):
    try:
        index_element=list.index(element)
        return index_element
    except ValueError:
        return -1
