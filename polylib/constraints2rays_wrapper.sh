#!/bin/bash

cat /tmp/.polylib-main | constraints2rays | sed 's~[ ][ ]*~ ~g' | sed 's~^ ~~' | sed 's~ $~~'
