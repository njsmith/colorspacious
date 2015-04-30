# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

class Trigger(object):
    def __init__(self):
        self._callbacks = set()

    def add_callback(self, f):
        self._callbacks.add(f)
        # Always call it immediately -- this is always legal (b/c as soon as
        # you call add_callback you have to be prepared for changes to
        # happen), saves having to explicitly call refresh methods in every
        # __init__, and flushes out bugs.
        f()

    def remove_callback(self, f):
        self._callbacks.remove(f)

    def fire(self):
        for f in self._callbacks:
            f()
