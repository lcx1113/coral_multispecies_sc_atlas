"""Lossless display compression of adjacent shared-prefix integer gene names."""
import re


def split_numeric_suffix(token):
    match = re.fullmatch(r'([A-Za-z][A-Za-z0-9]*?)([0-9]+)', token)
    return match.groups() if match else None


def simplify(annotation):
    head, separator, body = annotation.rpartition(':')
    if not separator:
        raise ValueError('Expected an annotated family prefix: ' + annotation)
    tokens = body.split('/')
    result, start = [], 0
    while start < len(tokens):
        parsed = split_numeric_suffix(tokens[start])
        end = start + 1
        if parsed:
            while end < len(tokens):
                following = split_numeric_suffix(tokens[end])
                if not following or following[0] != parsed[0]:
                    break
                end += 1
        result.append(tokens[start])
        result.extend(split_numeric_suffix(tokens[i])[1] for i in range(start + 1, end))
        start = end
    return head + separator + '/'.join(result)


def expand(annotation):
    head, separator, body = annotation.rpartition(':')
    result, active_stem = [], None
    for token in body.split('/'):
        if re.fullmatch(r'[0-9]+', token):
            if active_stem is None:
                raise ValueError('Number without a preceding full name: ' + annotation)
            result.append(active_stem + token)
        else:
            result.append(token)
            parsed = split_numeric_suffix(token)
            active_stem = parsed[0] if parsed else None
    return head + separator + '/'.join(result)


def truncate(annotation):
    return annotation[:65] + '...' if len(annotation) > 68 else annotation
