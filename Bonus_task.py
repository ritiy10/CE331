def new_series_map_number(lat, lon):

    # 4° latitude band
    band = int(lat // 4)
    first_letter = chr(ord('A') + band)

    # 6° longitude zone
    zone = int((lon + 180) // 6) + 1

    # Position within 4° × 6° block
    lon_start = (zone - 31) * 6
    col = int(lon - lon_start)

    lat_top = (band + 1) * 4
    row = int(lat_top - lat)

    # 1° × 1° sheet
    index = row * 6 + col
    sheet_letter = chr(ord('A') + index)

    # 15' × 15' sheet
    lat_start = int(lat)
    lon_start = int(lon)

    row15 = int((lat_start + 1 - lat) / 0.25)
    col15 = int((lon - lon_start) / 0.25)

    number = col15 * 4 + row15 + 1

    return f"{first_letter}{zone}{sheet_letter}{number}"


# India Gate
lat = 28.6129
lon = 77.2296

print(new_series_map_number(lat, lon))
