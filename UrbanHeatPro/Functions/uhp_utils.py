from typing import Any


def nested_get(input_dict: dict, nested_key: list, not_there="NotThere") -> Any:
    """
    Get value from nested dictionary given a list of nested keys and return default for missing keys.
    
    :param input_dict: Nested dictionary.
    :param nested_key: List of keys to access from outer dictionary to inner.
    :param not_there: Default value to be returned in case key does not exist
    :return: Value for the given list of keys. If a key is not available, then returns the default value.

    """
    internal_dict_value = input_dict
    for k in nested_key:
        internal_dict_value = internal_dict_value.get(k, not_there)
        if internal_dict_value is not_there:
            return not_there
    return internal_dict_value


def access_config_or_default(config: dict, default_config: dict, nested_key) -> Any:
    """
    Get the value for a given configuration parameter and use the default configuration value if the configuration
    parameter is not in the configuration there.
    
    :param config: Nested dictionary.
    :param default_config:
    :param nested_key: List of keys to access from outer dictionary to inner.
    :return: Configuration value for the given list of keys. If a key is not available, raises a KeyError.
    """
    lookup_res = nested_get(config, nested_key, nested_get(default_config, nested_key, "NotThere"))
    if lookup_res == "NotThere":
        raise KeyError("There was no default value in the config for", nested_key)
    return lookup_res
