# Photometric Recommender Pipeline (prep)

`prep` is a recommendation pipeline for new Type Ia Supernovae and was made to recommend potential triggers for the HST Dust project. It's data sources includes [ANTARES](https://antares.noirlab.edu/), [ALeRCE](https://alerce.science/), and [YSE](https://yse.ucsc.edu/). It has the ability to organize observations, perform as SALT3 fit and submit the information to a Slack channel or generate a csv file. This pipeline is currently in development.

## Adding `auth.py`

You will need an `auth.py` file to access [YSE-PZ](https://ziggy.ucolick.org/yse/dashboard/) and to post to Slack using their [API](https://api.slack.com/). It should have the path `./prep/auth.py` and should have the format,

```python
login = 'YOUR YSE-PZ LOGIN'
password = 'YOUR YSE-PZ PASSWORD'
toku = 'SLACK BOT/USER TOKEN'
```
You don't need the YSE-PZ credentials if you are just calling observations from ANTARES or ALeRCE. You can bypass this by setting `login` and `password` to `None` or `''`.
