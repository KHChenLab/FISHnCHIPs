from scipy import ndimage
import numpy as np

def _upsampled_dft(data, upsampled_region_size,
                   upsample_factor=1, axis_offsets=None):

    # if people pass in an integer, expand it to a list of equal-sized sections
    if not hasattr(upsampled_region_size, "__iter__"):
        upsampled_region_size = [upsampled_region_size, ] * data.ndim
    else:
        if len(upsampled_region_size) != data.ndim:
            raise ValueError("shape of upsampled region sizes must be equal "
                             "to input data's number of dimensions.")

    if axis_offsets is None:
        axis_offsets = [0, ] * data.ndim
    else:
        if len(axis_offsets) != data.ndim:
            raise ValueError("number of axis offsets must be equal to input "
                             "data's number of dimensions.")

    im2pi = 1j * 2 * np.pi

    dim_properties = list(zip(data.shape, upsampled_region_size, axis_offsets))

    for (n_items, ups_size, ax_offset) in dim_properties[::-1]:
        kernel = ((np.arange(ups_size) - ax_offset)[:, None]
                  * np.fft.fftfreq(n_items, upsample_factor))
        kernel = np.exp(-im2pi * kernel)

        # Equivalent to:
        #   data[i, j, k] = kernel[i, :] @ data[j, k].T
        data = np.tensordot(kernel, data, axes=(1, -1))
    return data

def register_translation(src_image: np.ndarray,
                         target_image: np.ndarray,
                         upsample_factor: int = 1,
                         space: str = "real",
                         return_error: bool = True,
                         ):

    # images must be the same shape
    if src_image.shape != target_image.shape:
        raise ValueError("Error: images must be same size for "
                         "register_translation")

    # assume complex data is already in Fourier space
    if space.lower() == 'fourier':
        src_freq = src_image
        target_freq = target_image
    # real data needs to be fft'd.
    elif space.lower() == 'real':
        src_freq = np.fft.fftn(src_image)
        target_freq = np.fft.fftn(target_image)
    else:
        raise ValueError("Error: register_translation only knows the \"real\" "
                         "and \"fourier\" values for the ``space`` argument.")

    # Whole-pixel shift - Compute cross-correlation by an IFFT
    shape = src_freq.shape
    image_product = src_freq * target_freq.conj()
    cross_correlation = np.fft.ifftn(image_product)

    # Locate maximum
    maxima = np.unravel_index(np.argmax(np.abs(cross_correlation)),
                              cross_correlation.shape)
    midpoints = np.array([np.fix(axis_size / 2) for axis_size in shape])

    shifts = np.array(maxima, dtype=np.float64)
    shifts[shifts > midpoints] -= np.array(shape)[shifts > midpoints]

    # calculate the pixel-resolution maximum cross-correlation value
    CCmax_pixel = cross_correlation[maxima]

    if upsample_factor == 1:
        if return_error:
            src_amp = np.sum(np.abs(src_freq) ** 2) / src_freq.size
            target_amp = np.sum(np.abs(target_freq) ** 2) / target_freq.size
            CCmax = CCmax_pixel # this is height of cc peak

    # If upsampling > 1, then refine estimate with matrix multiply DFT
    else:
        # Initial shift estimate in upsampled grid
        shifts = np.round(shifts * upsample_factor) / upsample_factor
        upsampled_region_size = np.ceil(upsample_factor * 1.5)

        # Center of output array at dftshift + 1
        dftshift = np.fix(upsampled_region_size / 2.0)
        upsample_factor = np.array(upsample_factor, dtype=np.float64)
        normalization = (src_freq.size * upsample_factor ** 2)

        # Matrix multiply DFT around the current shift estimate
        sample_region_offset = dftshift - shifts * upsample_factor
        cross_correlation = _upsampled_dft(image_product.conj(),
                                           upsampled_region_size,
                                           upsample_factor,
                                           sample_region_offset).conj()
        cross_correlation /= normalization

        # Locate maximum and map back to original pixel grid
        maxima = np.unravel_index(np.argmax(np.abs(cross_correlation)),
                                  cross_correlation.shape)
        CCmax = cross_correlation[maxima]

        maxima = np.array(maxima, dtype=np.float64) - dftshift

        shifts = shifts + maxima / upsample_factor

        if return_error:
            src_amp = _upsampled_dft(src_freq * src_freq.conj(),
                                     1, upsample_factor)[0, 0]
            src_amp /= normalization
            target_amp = _upsampled_dft(target_freq * target_freq.conj(),
                                        1, upsample_factor)[0, 0]
            target_amp /= normalization

    # For singleton dimensions, the shift calculated has no effect.
    # Set to zero.
    for dim in range(src_freq.ndim):
        if shape[dim] == 1:
            shifts[dim] = 0

    if return_error:
        # return (shifts,
        #         _compute_error(CCmax, src_amp, target_amp),
        #         _compute_phasediff(CCmax)
        #         )

        return (shifts,
                np.abs(CCmax),
                np.abs(CCmax_pixel)
                )
    else:
        return shifts

def register_slice(ref_slice, current_slice):
    ref_slice_fourier = np.fft.fftn(ref_slice)
    current_slice_fourier = np.fft.fftn(current_slice)
    (shifts, fine_error,  pixel_error) = register_translation(ref_slice_fourier,
                                                              current_slice_fourier,
                                                              upsample_factor=50,
                                                              space="fourier")
    registered_slice = np.fft.ifftn(ndimage.fourier_shift(current_slice_fourier, shifts))
    return registered_slice.real, shifts

def register_slice_with_shift(current_slice, shifts=(0,0)):

    current_slice_fourier = np.fft.fftn(current_slice)

    return np.fft.ifftn(ndimage.fourier_shift(current_slice_fourier, shifts)).real

def merge_image(image1, image2):

    # temporary RGB array to compare the reference and offset images
    rgb_temp = np.zeros(
        (image1.shape[0], image1.shape[1], 3),
        dtype=np.float32,
    )

    ref_temp = image1
    ref_min = np.min(ref_temp)
    ref_max = np.percentile(ref_temp, 99.98)
    ref_norm = (ref_temp - ref_min) / (ref_max - ref_min)
    rgb_temp[:, :, 0] = ref_norm

    reg_temp = image2
    reg_min = np.min(reg_temp)
    reg_max = np.percentile(reg_temp, 99.99)
    rgb_temp[:, :, 1] = (reg_temp - reg_min) / (reg_max - reg_min)

    return rgb_temp
