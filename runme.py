
from UrbanHeatPro import run_uhp

if __name__ == '__main__':

    # change the settings in settings/uhp_settings_custom.yaml
    # and then call run_uhp(settings_file='settings/uhp_settings_custom.yaml')

    # run_uhp(settings_file='settings/uhp_settings_example.yaml') for an example with the provided data for Unterhaching
    run_uhp(settings_file='settings/uhp_settings_example.yaml')
