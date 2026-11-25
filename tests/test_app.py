from streamlit.testing.v1 import AppTest
from pathlib import Path
APP = Path(__file__).resolve().parents[1] / 'app.py'


def test_app_default_and_changed_inputs():
    app=AppTest.from_file(APP,default_timeout=60).run()
    assert not app.exception
    assert not app.error
    assert any(m.value=='2.0303 м' for m in app.metric)
    app.selectbox(key='section').set_value('trapezoidal').run()
    assert not app.exception
    assert any('Длина' in w.value for w in app.warning)
    app.number_input(key='q').set_value(1.).run()
    assert not app.exception
    assert app.error


def test_scenarios_and_json_download_present():
    app=AppTest.from_file(APP,default_timeout=60).run()
    app.button[0].click().run()
    assert not app.exception
    assert len(app.session_state['batch_result'][1])==3


