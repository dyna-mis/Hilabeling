#ifndef ACLABELTESTER_H
#define ACLABELTESTER_H

#include <iostream>

#include "qgsproject.h"
#include "qgsdatasourceuri.h"
#include "qgsvectorlayer.h"
#include "qgsapplication.h"
#include "qgsvectorlayerlabelprovider.h"
#include "qgslabelfeature.h"
#include "qgsmaprenderersequentialjob.h"
#include "qgsfontutils.h"
#include "qgsvectorlayerlabeling.h"
#include "qgsvectorlayerlabelprovider.h"
#include "geometry/qgsrectangle.h"
#include "qgsauxiliarystorage.h"
#include <fstream>

#include <QMap>
#include <QDebug>
#include <vector>
#include "modification.h"
#include <iostream>
using namespace std;
namespace lt {
    json node;
    const QString END(" end"); 
    bool myinitial = true;
    unordered_map<int, int> my_solution_prev; 
    QTextStream& qStdOut()
    {
        static QTextStream ts( stdout );
        return ts;
    }
    // test name
    //(aaa...a) same algorithms 
    const QString SAME("same");
    //(ehhh...h) exact algorithms and then heuristic
    const QString CHANGE("change");
    //(ehhhh ehhhh,    ehhhh) exact algorithm before each period
    const QString INTER("inter");
    // runs
    int runs = 20;
    // algorithms
    typedef QPair<QString, QgsLabelingEngineSettings::Search> AlgorithmPair;
    const AlgorithmPair FALP = {"falp", QgsLabelingEngineSettings::Falp};
    const AlgorithmPair CHAIN = {"chain", QgsLabelingEngineSettings::Chain};
    const AlgorithmPair POPMUSIC = {"popmusic", QgsLabelingEngineSettings::Popmusic_Tabu_Chain};
    const AlgorithmPair SIMPLE = {"simple", QgsLabelingEngineSettings::SIMPLE};
    const AlgorithmPair MIS = {"mis", QgsLabelingEngineSettings::MIS};
    const AlgorithmPair KAMIS = {"kamis", QgsLabelingEngineSettings::KAMIS};
    const AlgorithmPair MAXHS = {"maxhs", QgsLabelingEngineSettings::MAXHS};
    QVector<AlgorithmPair> ALGORITHMS;
    QVector<AlgorithmPair> run_algorithms;
    // datasets
    QVector<QString> dataSets;
    const QString AUSTRIA("austria");
    const QString VIENNATRAM("vienna_subway_tram");
    const QString VIENNA("vienna_transport");
    class LabelTester
    {
        public:
              LabelTester(const QString& dataset, const QDir& base_path) :
           // LabelTester(const QString& dataset, const QDir& base_path = QDir("/home/fabian/clean_home/research/HILabeling/hilabeling_experiments/")) :
                ext_shapefile(".shp"),
                ext_image(".png"),
                map_size(5000, 5000),
                map_dpi(300),
                dataset(dataset),
                base_path(base_path),
                statistics_path(base_path.filePath("statistics")),
                shapefile_path(base_path.filePath("shapefiles")),
                image_path(base_path.filePath("images"))
            {
                Q_ASSERT(base_path.exists());
                Q_ASSERT(statistics_path.exists());
                Q_ASSERT(shapefile_path.exists());
                Q_ASSERT(image_path.exists());

                dataset_path = shapefile_path.filePath(dataset + ext_shapefile);
                Q_ASSERT(QFile(dataset_path).exists());
            }
            int data_test(const QString test)
            {
                my_solution_prev.clear();
                myinitial= true;
                run_test(test);
                return 0;
            }
            void prepare(){
                assert(myinitial);
                // Create layer & map
                qStdOut() << "Creating layer for dataset " << QFile(dataset_path).fileName() << endl;
                layer = new QgsVectorLayer(dataset_path, "", "ogr");
                QgsProject::instance()->addMapLayer(layer);
                Q_ASSERT(layer->isValid());
                qStdOut() << "Setting up map" << endl;
                setup_map();

                // Setting up PAL
                qStdOut() << "Setting up PAL" << endl;
                setup_pal();
                // Rendernig map once to save time later
                qStdOut() << "Rendering map" << endl;
                render_map();
            }
            bool run_test(const QString test)
            {
                myinitial = true;
               // Add field to track modifications
                // Test loop
                for(int i = 0; i < runs; ++i)
                {
//+++++++++++++++++++++++++++++++++MOdifucation+++++++++++++++++++++++++++++++++++++++++++++++++
                    if(i > 0){
                        myinitial = false;
                        Modification m;
                        vector<int> solutionIds;
                        m.modify(layer);
                        m.applyModification(fixExpressions, sizeExpressions);
                        //m.print();initial
                        // +++++++++++++++++apply Modification to layer +++++++++++++++++
                        if(debugTrue){
                            //vl->printLayer();
                        }
                        //------------------apply Modification to palSettings------------------
                       // qStdOut() <<"fixExpressions: " <<fixExpressions<< endl;
                       // qStdOut() <<"sizeExpressions: " <<sizeExpressions<< endl;
                        pal_settings.dataDefinedProperties().setProperty(QgsPalLayerSettings::Show,
                                    QgsProperty::fromExpression( fixExpressions + " end"));
                        pal_settings.dataDefinedProperties().setProperty(QgsPalLayerSettings::Size,
                                    QgsProperty::fromExpression( sizeExpressions + " end"));
                    }
                    else{
                        pal_settings.dataDefinedProperties().setProperty(QgsPalLayerSettings::Show,
                                    QgsProperty::fromExpression( ""));
                        pal_settings.dataDefinedProperties().setProperty(QgsPalLayerSettings::Size,
                                    QgsProperty::fromExpression("" ));
                           assert(myinitial);
                    }

//-------------------------------Modification-----------------------------------------------------                    

                    // Setting up layer & creating provider
                    layer->setLabelsEnabled(true);
                    QgsVectorLayerLabelProvider *provider = new QgsVectorLayerLabelProvider(layer, "test", true, &pal_settings);
                    layer->setLabeling(new QgsVectorLayerSimpleLabeling(pal_settings));

                    AlgorithmPair algorithm = run_algorithms[i];
                    qStdOut() << "Running test with algorithm " << algorithm.first << endl;

                    //Get rendering context
                    QImage labeled_image = map_image.copy();
                    QPainter p(&labeled_image);
                    QgsRenderContext context = QgsRenderContext::fromMapSettings(map_settings);
                    context.setPainter(&p);

                    // Configuring algorithm
                    qStdOut() << "Setting algorithm" << endl;
                    QgsLabelingEngineSettings engine_settings;
                    engine_settings.setFlag(QgsLabelingEngineSettings::UsePartialCandidates, false);
                    engine_settings.setSearchMethod(algorithm.second);
                    map_settings.setLabelingEngineSettings(engine_settings);

                    // Create Labeling Engine
                    qStdOut() << "Creating labeling engine" << endl;
                    QgsLabelingEngine engine;
                    engine.setMapSettings(map_settings);
                    engine.addProvider(provider);

                    // Create one labeling with algorithm
                    create_labeling_with_algorithm(engine, context, map_settings,i);

                    p.end();

                    // Store image
                    QString save_to_image_path(image_path.filePath(test+ "_"+dataset +"_" + QString::number(i) +  ext_image));
                    qStdOut() << "Saving labeled map to: " << save_to_image_path << endl;
                    labeled_image.save(save_to_image_path);
                }
                QString save_to_statistics_path(statistics_path.filePath(test+ "_"+dataset +".json"));
                std::fstream fs;
                fs.open (save_to_statistics_path.toUtf8().constData(),std::fstream::out);
                if (fs.is_open()){
                        fs<< node;
                        fs.close();
                }
                else{
                    std::cout << "Error opening file";
                }
                sizeExpressions =" case ";
                fixExpressions = " case ";
                node.clear();
                return true;
            }
        private:
            /*
             * One labeling with one algorithm for the given map and provider is created
             */
            QgsLabelingResults* create_labeling_with_algorithm(QgsLabelingEngine& engine, QgsRenderContext& context,
                                                QgsMapSettings& map_settings,int i)
            {
                // Run engine
                qStdOut() << "Running labeling engine..." << endl;
                test::Performance performance;
                engine.run(context,performance,myinitial,my_solution_prev);
                if(i ==0) assert(performance.remainingLabels == 0);
                node["Round " + std::to_string(i)] = performance.convertJSON();
                qStdOut() << "Labeled " << engine.results()->labelsWithinRect(map_settings.layers().first()->extent()).size() << " features" << endl;
                return engine.results();
            }

            bool setup_map()
            {
                map_settings.setOutputSize(map_size);
                map_settings.setExtent(layer->extent());
                map_settings.setLayers(QList<QgsMapLayer *>() << layer);
                map_settings.setOutputDpi(map_dpi);

                return true;
            }

            bool setup_pal()
            {
                pal_settings.fieldName = QStringLiteral("name");

                QgsTextFormat format;
                format.setFont(QgsFontUtils::getStandardTestFont(QStringLiteral("Bold")).family());
                format.setSize(10);
                format.setNamedStyle(QStringLiteral("Bold"));
                format.setColor(QColor(200, 0, 200));
                pal_settings.setFormat(format);
                pal_settings.placement = QgsPalLayerSettings::OrderedPositionsAroundPoint;

                return true;
            }

            bool render_map()
            {
                // Get rendering context
                qStdOut() << "Get rendering context" << endl;
                QgsMapRendererSequentialJob job(map_settings);
                              cout<< "see what happens A";
                myinitial = true;
                              cout<< "see what happens B";
                job.start();
                cout<< "see what happens C";
                job.waitForFinished();
                map_image = job.renderedImage();

                return true;
            }
            QgsVectorLayer *layer;
            QgsMapSettings map_settings;
            QgsPalLayerSettings pal_settings;
            QImage map_image;
            QString sizeExpressions =" case ";
            QString fixExpressions = " case ";
            //QString sizeExpressions = QStringLiteral("case when \"name\"='Oberlaa' then 30");
            //QString fixExpressions = QStringLiteral( "case when \"name\"='Leopoldau' then '#eeeeee'" );
            QString ext_shapefile;
            QString ext_image;
            QSize map_size;
            qint16 map_dpi;

            QString dataset;
            QString dataset_path;

            QDir base_path;
            QDir statistics_path;
            QDir shapefile_path;
            QDir image_path;
    };

    void sameAlgoGen(AlgorithmPair initA){
        run_algorithms.clear();
        for(int i = 0; i < runs; i++){
            run_algorithms.push_back(initA);
        }
    }
    void changeAlgoGen(AlgorithmPair initA, AlgorithmPair updateA){
        run_algorithms.clear();
        run_algorithms.push_back(initA);
        for(int i = 0; i < runs-1; i++){
            run_algorithms.push_back(updateA);
        }
    }
    void interAlgoGen(AlgorithmPair interA, AlgorithmPair updateA){
        run_algorithms.clear();
        for(int i = 0; i < 4; i++){
            run_algorithms.push_back(interA);
            for(int j = 0; j < 4; j++){
                run_algorithms.push_back(updateA);
            }
        }
    }
    void prepare(){
        // prepare the datas;
        //dataSets.push_back(AUSTRIA);
        //dataSets.push_back(VIENNA);
        dataSets.push_back(VIENNATRAM);
        // prepare the  7 algorithms
        ALGORITHMS.push_back(FALP);
        ALGORITHMS.push_back(CHAIN);
        ALGORITHMS.push_back(POPMUSIC);
        ALGORITHMS.push_back(SIMPLE);
        ALGORITHMS.push_back(KAMIS);
        ALGORITHMS.push_back(MIS);
        ALGORITHMS.push_back(MAXHS);
                //qDebug()<< dataSets;
                //qDebug()<< ALGORITHMS;
    }
    int test(const QString& dataset){
        prepare();
        const QDir& base_path = QDir("/home/guangping/dev/cpp");
        qDebug()<< dataSets;
        for(auto data: dataSets){
                LabelTester tester(dataset, base_path);
                tester.prepare();
                // TEST A (aaa... aaa) <runs> times (default 20) with  same algorithm
                for(auto initA: ALGORITHMS){
                    sameAlgoGen(initA);
                    QString comment = SAME + "_"+initA.first;
                    tester.data_test(comment);
                }
                // TEST B (abbb... bbb)
                for(auto initA: ALGORITHMS){
                    for(auto updateA: ALGORITHMS){
                        if(initA.first == updateA.first) continue;
                        changeAlgoGen(initA, updateA);
                        QString comment = CHANGE + "_"+initA.first+"_"+updateA.first;
                        tester.data_test(comment);
                    }
                }
                // TEST C (abbbb... abbbb)
                for(auto initA: ALGORITHMS){
                    for(auto updateA: ALGORITHMS){
                        if(initA.first == updateA.first) continue;
                        interAlgoGen(initA, updateA);
                        QString comment = INTER + "_"+initA.first+"_"+updateA.first;
                        tester.data_test(comment);
                    }
                }
        }
        return 0;
    }
}


#endif // ACLABELTESTER_H