
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;
import java.util.Timer;
import java.util.TimerTask;


public class ImageDisplay {

    JFrame frame;
    JLabel lbIm1;
    JLabel lbOriginal;
    BufferedImage imgOne;

    // Modify the height and width values here to read and display an image with
    // different dimensions.
    int width = 512;
    int height = 512;

    private int[][][][] colorBlocks;
    private int currentBlockX = 0;
    private int currentBlockY = 0;
    private BufferedImage displayImage;

    //Input parameters
    int quantizationLevel;
    String deliveryMode;
    int latency;

    //Mode 2
    private int currentCoefficientIndex = 0;
    private int currentBitDepth = 1;

    //Precompute values
    private double[][] cosValues = new double[8][8];


    private void precomputeCosValues() {
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                cosValues[i][j] = Math.cos((2 * i + 1) * j * Math.PI / 16);

            }
        }
    }

    /**
     * Read Image RGB
     * Reads the image of given width and height at the given imgPath into the provided BufferedImage.
     */
    private void readImageRGB(int width, int height, String imgPath, BufferedImage img) {
        try {
            int frameLength = width * height * 3;

            File file = new File(imgPath);
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            raf.seek(0);

            long len = frameLength;
            byte[] bytes = new byte[(int) len];

            raf.read(bytes);

            int ind = 0;
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    byte a = 0;
                    byte r = bytes[ind];
                    byte g = bytes[ind + height * width];
                    byte b = bytes[ind + height * width * 2];

                    int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
                    //int pix = ((a << 24) + (r << 16) + (g << 8) + b);
                    img.setRGB(x, y, pix);
                    ind++;
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void showIms(String[] args) {

        // Read in the specified image
        imgOne = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        displayImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        readImageRGB(width, height, args[0], imgOne);

        // Use label to display the image
        frame = new JFrame();
        GridLayout gLayout = new GridLayout(1, 2);
        frame.getContentPane().setLayout(gLayout);

        lbIm1 = new JLabel(new ImageIcon(displayImage));
        lbOriginal = new JLabel(new ImageIcon(imgOne));

        GridBagConstraints c = new GridBagConstraints();
        c.anchor = GridBagConstraints.CENTER;
        c.weightx = 0.5;

        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridx = 0;
        c.gridy = 1;

        frame.getContentPane().add(lbOriginal, c);
        frame.getContentPane().add(lbIm1, c);

        frame.pack();
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Input parameters
        quantizationLevel = Integer.parseInt(args[1]);
        deliveryMode = args[2];
        latency = Integer.parseInt(args[3]);
        if (latency == 0) latency++; //Giving 1 millisecond when latency is 0 as the timer expects positive integer

        //Precompute
        precomputeCosValues();

        //Encode
        encodeImage(imgOne, quantizationLevel);

        //Decode and display
        decodeAndDisplayBasedOnDeliveryMode();
    }

    private void decodeAndDisplayBasedOnDeliveryMode() {
        switch (deliveryMode) {
            case "1":
                displaySequentialBlocks();
                break;
            case "2":
                displayProgressiveBlocks();
                break;
            case "3":
                displayProgressiveBitWiseBlocks();
                break;
            default:
                System.out.println("Invalid delivery mode selected.");
                break;
        }
    }

    // New method to initialize and start the block display sequence
    private void displaySequentialBlocks() {
        Timer timer = new Timer();

        timer.scheduleAtFixedRate(new TimerTask() {
            @Override
            public void run() {
                if (currentBlockY < colorBlocks[0].length) {
                    if (currentBlockX < colorBlocks[0][currentBlockY].length) {
                        // Extract the current block for each color channel
                        int[][] blockR = new int[8][8];
                        int[][] blockG = new int[8][8];
                        int[][] blockB = new int[8][8];
                        for (int i = 0; i < 64; i++) {
                            blockR[i / 8][i % 8] = colorBlocks[0][currentBlockY][currentBlockX][i];
                            blockG[i / 8][i % 8] = colorBlocks[1][currentBlockY][currentBlockX][i];
                            blockB[i / 8][i % 8] = colorBlocks[2][currentBlockY][currentBlockX][i];
                        }

                        processAndDisplayBlocks(blockR, blockG, blockB);

                        //Display block-wise
                        lbIm1.setIcon(new ImageIcon(displayImage));
                        frame.repaint();
                        currentBlockX++;
                    } else {
                        currentBlockX = 0;
                        currentBlockY++;
                    }
                } else {
                    // Cancel the timer when all blocks have been displayed
                    timer.cancel();
                }
            }
        }, 0, latency);
    }

    private void processAndDisplayBlocks(int[][] blockR, int[][] blockG, int[][] blockB) {
        //Dequantize blocks
        dequantizeBlock(blockR);
        dequantizeBlock(blockG);
        dequantizeBlock(blockB);

        // Apply IDCT
        applyIDCTToBlock(blockR);
        applyIDCTToBlock(blockG);
        applyIDCTToBlock(blockB);

        // Draw the processed blocks onto the display image
        drawBlock(blockR, blockG, blockB, currentBlockX, currentBlockY);
    }


    private void displayProgressiveBlocks() {
        Timer timer = new Timer();

        timer.scheduleAtFixedRate(new TimerTask() {
            @Override
            public void run() {
                // Decode and display using increasing number of coefficients
                for (int by = 0; by < colorBlocks[0].length; by++) {
                    for (int bx = 0; bx < colorBlocks[0][by].length; bx++) {
                        int[][] blockR = extractBlock(colorBlocks[0], bx, by);
                        int[][] blockG = extractBlock(colorBlocks[1], bx, by);
                        int[][] blockB = extractBlock(colorBlocks[2], bx, by);
                        processAndDisplayBlockCoeffWise(blockR, blockG, blockB, bx, by, currentCoefficientIndex);
                    }
                }

                lbIm1.setIcon(new ImageIcon(displayImage));
                frame.repaint();

                currentCoefficientIndex++;
                if (currentCoefficientIndex > 63) { // Finished all coefficients
                    timer.cancel();
                }
            }
        }, 0, latency);
    }

    private int[][] extractBlock(int[][][] blocks, int bx, int by) {
        int[][] block = new int[8][8];
        for (int i = 0; i < 64; i++) {
            block[i / 8][i % 8] = blocks[by][bx][i];
        }
        return block;
    }

    private void processAndDisplayBlockCoeffWise(int[][] blockR, int[][] blockG, int[][] blockB,
                                                 int blockX, int blockY, int coeffIndex) {
        //Dequantize first
        dequantizeBlock(blockR);
        dequantizeBlock(blockG);
        dequantizeBlock(blockB);

        // Set coefficients after the current index to zero
        zeroOutCoefficientsAfterIndex(blockR, coeffIndex);
        zeroOutCoefficientsAfterIndex(blockG, coeffIndex);
        zeroOutCoefficientsAfterIndex(blockB, coeffIndex);

        // Apply IDCT
        applyIDCTToBlock(blockR);
        applyIDCTToBlock(blockG);
        applyIDCTToBlock(blockB);

        // Draw the processed blocks onto the display image
        drawBlock(blockR, blockG, blockB, blockX, blockY);
    }

    private void zeroOutCoefficientsAfterIndex(int[][] block, int index) {
        // ZigZag Order for a 8x8 Block
        int[][] zigZagOrder = {
                {0, 1, 5, 6, 14, 15, 27, 28},
                {2, 4, 7, 13, 16, 26, 29, 42},
                {3, 8, 12, 17, 25, 30, 41, 43},
                {9, 11, 18, 24, 31, 40, 44, 53},
                {10, 19, 23, 32, 39, 45, 52, 54},
                {20, 22, 33, 38, 46, 51, 55, 60},
                {21, 34, 37, 47, 50, 56, 59, 61},
                {35, 36, 48, 49, 57, 58, 62, 63}
        };

        for (int i = index + 1; i < 64; i++) {
            // Convert the ZigZag index to x and y coordinates
            int x = -1, y = -1;
            search:
            for (int row = 0; row < zigZagOrder.length; row++) {
                for (int col = 0; col < zigZagOrder[row].length; col++) {
                    if (zigZagOrder[row][col] == i) {
                        x = row;
                        y = col;
                        break search;
                    }
                }
            }

            // If the coordinate is found, zero out the coefficient
            if (x != -1 && y != -1) {
                block[x][y] = 0;
            }
        }
    }


    private void displayProgressiveBitWiseBlocks() {
        Timer timer = new Timer();
        timer.scheduleAtFixedRate(new TimerTask() {
            @Override
            public void run() {
                // Iterate over all blocks
                for (int by = 0; by < colorBlocks[0].length; by++) {
                    for (int bx = 0; bx < colorBlocks[0][by].length; bx++) {
                        int[][] blockR = extractBlock(colorBlocks[0], bx, by);
                        int[][] blockG = extractBlock(colorBlocks[1], bx, by);
                        int[][] blockB = extractBlock(colorBlocks[2], bx, by);

                        // Process the blocks with the current bit depth
                        processAndDisplayBlockBitwise(blockR, blockG, blockB, bx, by, currentBitDepth);
                    }
                }

                lbIm1.setIcon(new ImageIcon(displayImage));
                frame.repaint();

                currentBitDepth++;
                if (currentBitDepth > 32) {
                    timer.cancel();
                }
            }
        }, 0, latency);
    }

    private void processAndDisplayBlockBitwise(int[][] blockR, int[][] blockG, int[][] blockB,
                                               int blockX, int blockY, int bitDepth) {
        // First, dequantize
        dequantizeBlock(blockR);
        dequantizeBlock(blockG);
        dequantizeBlock(blockB);

        // Modify each block for the current bit depth
        modifyBlockForBitDepth(blockR, bitDepth);
        modifyBlockForBitDepth(blockG, bitDepth);
        modifyBlockForBitDepth(blockB, bitDepth);

        // Apply IDCT
        applyIDCTToBlock(blockR);
        applyIDCTToBlock(blockG);
        applyIDCTToBlock(blockB);

        // Draw the block
        drawBlock(blockR, blockG, blockB, blockX, blockY);
    }

    private void modifyBlockForBitDepth(int[][] block, int bitDepth) {
        for (int i = 0; i < block.length; i++) {
            for (int j = 0; j < block[i].length; j++) {
                block[i][j] = approximateCoefficient(block[i][j], bitDepth);
            }
        }
    }

    private int approximateCoefficient(int value, int bitDepth) {
        boolean isNegative = value < 0;
        if (isNegative) {
            value = -value;
        }
        int rebuiltValue = 0;
        int bitValue = 1 << 31;
        for (int i = 0; i < bitDepth; i++) {
            if ((value & bitValue) != 0) {
                rebuiltValue |= bitValue;
            }
            bitValue >>>= 1;
        }
        if (isNegative) {
            rebuiltValue = -rebuiltValue;
        }
        return rebuiltValue;
    }

    //Method do draw a given decoded block
    private void drawBlock(int[][] blockR, int[][] blockG, int[][] blockB, int blockX, int blockY) {
        for (int y = 0; y < 8; y++) {
            for (int x = 0; x < 8; x++) {
                // Ensure the pixel value is within the 0-255 range after processing
                int red = blockR[y][x];
                int green = blockG[y][x];
                int blue = blockB[y][x];

                // Combine the RGB values into a single integer
                int rgb = (red << 16) | (green << 8) | blue;
                displayImage.setRGB(blockX * 8 + x, blockY * 8 + y, rgb);
            }
        }
    }

    // Method to encode image by applying DCT and Quantization
    private void encodeImage(BufferedImage img, int quantizationLevel) {
        int blocksX = img.getWidth() / 8;
        int blocksY = img.getHeight() / 8;
        colorBlocks = new int[3][blocksY][blocksX][64];

        for (int by = 0; by < blocksY; by++) {
            for (int bx = 0; bx < blocksX; bx++) {
                for (int y = 0; y < 8; y++) {
                    for (int x = 0; x < 8; x++) {
                        int pixel = img.getRGB(bx * 8 + x, by * 8 + y);
                        colorBlocks[0][by][bx][y * 8 + x] = ((pixel >> 16) & 0xff);
                        colorBlocks[1][by][bx][y * 8 + x] = ((pixel >> 8) & 0xff);
                        colorBlocks[2][by][bx][y * 8 + x] = (pixel & 0xff);
                    }
                }

                // Apply DCT and Quantization to each block for each color channel
                for (int i = 0; i < 3; i++) {
                    int[][] block = new int[8][8];
                    for (int j = 0; j < 64; j++) {
                        block[j / 8][j % 8] = colorBlocks[i][by][bx][j];
                    }

                    applyDCTToBlock(block);
                    quantizeBlock(block, quantizationLevel);

                    // Store the quantized block back into the colorBlocks array
                    for (int j = 0; j < 64; j++) {
                        colorBlocks[i][by][bx][j] = block[j / 8][j % 8];
                    }
                }
            }
        }
    }

    // DCT transformation for a block
    private void applyDCTToBlock(int[][] block) {
        int size = 8; // Block size for DCT
        double[][] temp = new double[size][size];

        for (int u = 0; u < size; u++) {
            for (int v = 0; v < size; v++) {
                double sum = 0.0;
                for (int x = 0; x < size; x++) {
                    for (int y = 0; y < size; y++) {
                        sum += block[x][y] *
                                cosValues[x][u] *
                                cosValues[y][v];
                    }
                }
                double alphaU = (u == 0) ? 1.0 / Math.sqrt(2) : 1.0;
                double alphaV = (v == 0) ? 1.0 / Math.sqrt(2) : 1.0;
                temp[u][v] = 0.25 * alphaU * alphaV * sum;
            }
        }

        for (int u = 0; u < size; u++) {
            for (int v = 0; v < size; v++) {
                block[u][v] = (int)temp[u][v];
            }
        }
    }

    // Method to quantize a block using a uniform quantization table
    private void quantizeBlock(int[][] block, int quantizationLevel) {
        int quantizationFactor = (int) Math.pow(2, quantizationLevel);
        for (int u = 0; u < block.length; u++) {
            for (int v = 0; v < block[u].length; v++) {
                block[u][v] = Math.round(block[u][v] / quantizationFactor);
            }
        }
    }

    // Method to dequantize a block using a uniform quantization level
    private void dequantizeBlock(int[][] block) {
        int quantizationFactor = (int) Math.pow(2, quantizationLevel);
        for (int u = 0; u < block.length; u++) {
            for (int v = 0; v < block[u].length; v++) {
                block[u][v] *= quantizationFactor;
            }
        }
    }

    // Method to apply IDCT to a block (8x8)
    private void applyIDCTToBlock(int[][] block) {
        int size = 8;
        double[][] temp = new double[size][size];
        double cu, cv, sum;

        for (int x = 0; x < size; x++) {
            for (int y = 0; y < size; y++) {
                sum = 0.0;
                for (int u = 0; u < size; u++) {
                    for (int v = 0; v < size; v++) {
                        cu = (u == 0) ? 1 / Math.sqrt(2) : 1.0;
                        cv = (v == 0) ? 1 / Math.sqrt(2) : 1.0;
                        sum += cu * cv * block[u][v] *
                                cosValues[x][u] *
                                cosValues[y][v];
                    }
                }
                temp[x][y] = sum / 4;
            }
        }
        for (int x = 0; x < size; x++) {
            for (int y = 0; y < size; y++) {
                block[x][y] = (int) Math.min(Math.max(temp[x][y], 0), 255);
            }
        }
    }

    public static void main(String[] args) {
        ImageDisplay ren = new ImageDisplay();
        ren.showIms(args);
    }

}
