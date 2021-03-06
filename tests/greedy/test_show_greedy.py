# -*- coding: utf-8 -*-
"""Test map visualization.

These tests don't check for anything; we're just ensuring
the code runs without raising errors.

"""
import matplotlib

matplotlib.use("Agg")

import starry
import numpy as np
import os


def test_show(mp4=False):
    map = starry.Map(ydeg=1, udeg=1)
    map.show(file="tmp.pdf", projection="ortho")
    os.remove("tmp.pdf")
    map.show(file="tmp.pdf", projection="rect")
    os.remove("tmp.pdf")
    if mp4:
        map.show(theta=np.linspace(0, 360, 10), file="tmp.mp4")
        os.remove("tmp.mp4")


def test_show_reflected(mp4=False):
    map = starry.Map(ydeg=1, udeg=1, reflected=True)
    map.show(file="tmp.pdf", projection="ortho")
    os.remove("tmp.pdf")
    map.show(file="tmp.pdf", projection="rect")
    os.remove("tmp.pdf")
    if mp4:
        map.show(theta=np.linspace(0, 360, 10), file="tmp.mp4")
        os.remove("tmp.mp4")


def test_show_rv(mp4=False):
    map = starry.Map(ydeg=1, udeg=1, rv=True)
    map.show(rv=True, file="tmp.pdf", projection="ortho")
    os.remove("tmp.pdf")
    map.show(rv=True, file="tmp.pdf", projection="rect")
    os.remove("tmp.pdf")
    if mp4:
        map.show(rv=True, theta=np.linspace(0, 360, 10), file="tmp.mp4")
        os.remove("tmp.mp4")


def test_show_ld():
    map = starry.Map(udeg=2)
    map.show(file="tmp.pdf")
    os.remove("tmp.pdf")


def test_system_show(mp4=False):
    pri = starry.Primary(starry.Map())
    sec = starry.Secondary(starry.Map(), porb=1.0)
    sys = starry.System(pri, sec)
    sys.show(0.1, file="tmp.pdf")
    os.remove("tmp.pdf")
    if mp4:
        sys.show([0.1, 0.2], file="tmp.mp4")
        os.remove("tmp.mp4")
    sys.show([0.1, 0.2], file="tmp.gif")
    os.remove("tmp.gif")


def test_system_rv_show(mp4=False):
    pri = starry.Primary(starry.Map(rv=True))
    sec = starry.Secondary(starry.Map(rv=True), porb=1.0)
    sys = starry.System(pri, sec)
    sys.show(0.1, file="tmp.pdf")
    os.remove("tmp.pdf")
    if mp4:
        sys.show([0.1, 0.2], file="tmp.mp4")
        os.remove("tmp.mp4")
